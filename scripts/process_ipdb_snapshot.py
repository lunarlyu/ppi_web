#!/usr/bin/env python3
"""Build compact browser data from an upstream IPDB snapshot."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

DEFAULT_SNAPSHOT = Path(
    "~/Downloads/Induced Proximity Database Project copy/public/ipdb-snapshot.json"
)
DEFAULT_TARGET_SCORES = Path(
    "~/Downloads/Induced Proximity Database Project copy/exports/riptac_target_scores.csv"
)
DEFAULT_LITERATURE_POCKETS = Path(
    "~/Downloads/Induced Proximity Database Project copy/src/ipdb/data/literature_pockets.tsv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build compact RIPTAC, protein-card, and ligand indexes.",
    )
    parser.add_argument(
        "snapshot",
        type=Path,
        nargs="?",
        default=DEFAULT_SNAPSHOT,
        help="Path to the upstream IPDB ipdb-snapshot.json file.",
    )
    parser.add_argument(
        "--ipdb-output",
        type=Path,
        default=Path("data/ipdb_processed.json"),
        help="Processed IPDB JSON consumed by the browser app.",
    )
    parser.add_argument(
        "--ligand-output",
        type=Path,
        default=Path("data/known_ligands.json"),
        help="Compact BindingDB known-ligand JSON consumed by the browser app.",
    )
    parser.add_argument(
        "--target-score-input",
        type=Path,
        default=DEFAULT_TARGET_SCORES,
        help="RIPTAC target score CSV with protein_id, cancer_type, and score.",
    )
    parser.add_argument(
        "--literature-pocket-input",
        type=Path,
        default=DEFAULT_LITERATURE_POCKETS,
        help="Literature pocket TSV with gene_symbol, pocket_type, and evidence_label.",
    )
    return parser.parse_args()


def display_path(path: Path) -> str:
    try:
        return f"~/{path.relative_to(Path.home())}"
    except ValueError:
        return str(path)


def clean_binding(binding: dict[str, Any]) -> dict[str, Any]:
    return {
        "assay_type": binding.get("assay_type"),
        "protac_role": binding.get("protac_role"),
        "source": binding.get("source"),
        "value_nm": binding.get("value_nm"),
    }


def clean_pocket(pocket: dict[str, Any]) -> dict[str, Any]:
    return prune(
        {
            "best_affinity_nm": pocket.get("best_affinity_nm"),
            "count": pocket.get("count"),
            "ligandability_score": pocket.get("ligandability_score"),
            "pocket_type": pocket.get("pocket_type"),
            "sources": pocket.get("sources"),
        },
    )


def clean_evidence(record: dict[str, Any]) -> dict[str, Any]:
    cleaned = {
        "best_binding": clean_binding(record["best_binding"])
        if isinstance(record.get("best_binding"), dict)
        else None,
        "compartments": record.get("compartments"),
        "dependency": record.get("dependency"),
        "driver": record.get("driver"),
        "expression": record.get("expression"),
        "pockets": [
            clean_pocket(pocket)
            for pocket in record.get("pockets", [])
            if isinstance(pocket, dict)
        ],
    }

    return prune(cleaned)


def clean_pair_row(row: dict[str, Any]) -> dict[str, Any]:
    return prune(
        {
            "assets": row.get("assets"),
            "cancer_type": row.get("cancer_type"),
            "left": normalize_symbol(row.get("left")),
            "metrics": row.get("metrics"),
            "modality": row.get("modality"),
            "right": normalize_symbol(row.get("right")),
            "score": row.get("score"),
        },
    )


def clean_target_score_row(row: dict[str, Any]) -> dict[str, Any]:
    return prune(
        {
            "cancer_type": row.get("cancer_type"),
            "score": parse_float(row.get("score")),
        },
    )


def clean_literature_pocket_row(row: dict[str, Any]) -> dict[str, Any]:
    return prune(
        {
            "evidence_label": row.get("evidence_label"),
            "ligandability_score": parse_float(row.get("ligandability_score")),
            "pocket_type": row.get("pocket_type"),
        },
    )


def load_target_scores(
    target_score_path: Path,
    source_proteins: dict[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    if not target_score_path.exists():
        return {}

    symbols_by_protein_id = {
        str(record.get("protein_id") or "").strip(): normalize_symbol(symbol)
        for symbol, record in source_proteins.items()
        if isinstance(record, dict) and record.get("protein_id")
    }
    rows_by_symbol: dict[str, list[dict[str, Any]]] = {}

    with target_score_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            protein_id = str(row.get("protein_id") or "").strip()
            symbol = symbols_by_protein_id.get(protein_id)
            cleaned = clean_target_score_row(row)
            if not symbol or "cancer_type" not in cleaned or "score" not in cleaned:
                continue

            rows_by_symbol.setdefault(symbol, []).append(cleaned)

    return {
        symbol: sorted(rows, key=lambda row: row["score"], reverse=True)[:3]
        for symbol, rows in rows_by_symbol.items()
    }


def load_literature_pockets(
    literature_pocket_path: Path,
    source_proteins: dict[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    if not literature_pocket_path.exists():
        return {}

    symbols_by_alias = protein_aliases(source_proteins)
    rows_by_symbol: dict[str, list[dict[str, Any]]] = {}

    with literature_pocket_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            symbol = symbols_by_alias.get(normalize_symbol(row.get("gene_symbol")))
            cleaned = clean_literature_pocket_row(row)
            if not symbol or "pocket_type" not in cleaned:
                continue

            rows_by_symbol.setdefault(symbol, []).append(cleaned)

    return {
        symbol: sorted(
            rows,
            key=lambda row: (
                -(row.get("ligandability_score") or 0),
                row.get("evidence_label") or "",
            ),
        )
        for symbol, rows in rows_by_symbol.items()
    }


def protein_aliases(source_proteins: dict[str, Any]) -> dict[str, str]:
    aliases: dict[str, str] = {}

    for raw_symbol, record in source_proteins.items():
        if not isinstance(record, dict):
            continue

        symbol = normalize_symbol(raw_symbol)
        if not symbol:
            continue

        for alias in (
            raw_symbol,
            record.get("hgnc_symbol"),
            record.get("protein_id"),
            *(str(record.get("uniprot_accession") or "").split("|")),
        ):
            normalized = normalize_symbol(alias)
            if normalized:
                aliases[normalized] = symbol

    return aliases


def parse_float(value: Any) -> float | None:
    if value in ("", None):
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def normalize_symbol(value: Any) -> str:
    return str(value or "").strip().upper()


def normalize_pair_key(left: str, right: str) -> str:
    return "::".join(sorted([normalize_symbol(left), normalize_symbol(right)]))


def prune(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: pruned
            for key, item in value.items()
            if (pruned := prune(item)) not in (None, [], {})
        }

    if isinstance(value, list):
        return [
            pruned
            for item in value
            if (pruned := prune(item)) not in (None, [], {})
        ]

    return value


def main() -> None:
    args = parse_args()
    snapshot_path = args.snapshot.expanduser()

    with snapshot_path.open("r", encoding="utf-8") as handle:
        snapshot = json.load(handle)

    source_proteins = snapshot.get("proteins", {})
    proteins = {}
    evidence = {}
    ligands = {}
    for raw_symbol, record in snapshot.get("evidence", {}).items():
        if not isinstance(record, dict):
            continue

        symbol = normalize_symbol(raw_symbol)
        binding = record.get("best_binding")
        if symbol and isinstance(binding, dict):
            ligands[symbol] = clean_binding(binding)

        cleaned_evidence = clean_evidence(record)
        if symbol and cleaned_evidence:
            evidence[symbol] = cleaned_evidence

    pair_rows = {}
    for key, pair in snapshot.get("pairs", {}).items():
        if not isinstance(pair, dict):
            continue

        protac = [
            clean_pair_row(row)
            for row in pair.get("protac", [])
            if isinstance(row, dict)
        ]
        riptac = [
            clean_pair_row(row)
            for row in pair.get("riptac", [])
            if isinstance(row, dict)
        ]
        compact_pair = prune({"protac": protac, "riptac": riptac})
        if compact_pair:
            left, _, right = str(key).partition("::")
            pair_rows[normalize_pair_key(left, right)] = compact_pair
            for row in protac + riptac:
                evidence.setdefault(row["left"], {})
                evidence.setdefault(row["right"], {})

    target_scores = load_target_scores(
        args.target_score_input.expanduser(),
        source_proteins,
    )
    for symbol, scores in target_scores.items():
        evidence.setdefault(symbol, {})["target_scores"] = scores

    literature_pockets = load_literature_pockets(
        args.literature_pocket_input.expanduser(),
        source_proteins,
    )
    for symbol, pockets in literature_pockets.items():
        evidence.setdefault(symbol, {})["literature_pockets"] = pockets

    for symbol in evidence:
        record = source_proteins.get(symbol, {})
        if not isinstance(record, dict):
            continue

        name = record.get("name")
        if name:
            proteins[symbol] = {"name": name}

    generated_at = datetime.now(UTC).isoformat()
    source_meta = {
        "generated_at": generated_at,
        "source_file": display_path(snapshot_path),
        "snapshot_generated_at": snapshot.get("generated_at"),
    }

    ipdb_output = {
        "evidence": dict(sorted(evidence.items())),
        "meta": {
            "description": (
                "Compact IPDB pairs, protein names, and evidence records "
                "for the STRING PPI star map."
            ),
            "evidence_count": len(evidence),
            "literature_pocket_count": sum(
                len(pockets) for pockets in literature_pockets.values()
            ),
            "pair_count": len(pair_rows),
            "protein_count": len(proteins),
            "target_score_count": len(target_scores),
            **source_meta,
        },
        "pairs": dict(sorted(pair_rows.items())),
        "proteins": dict(sorted(proteins.items())),
    }

    ligand_output = {
        "ligands": dict(sorted(ligands.items())),
        "meta": {
            "description": (
                "Gene-level known-ligand index derived from IPDB "
                "evidence.best_binding BindingDB records."
            ),
            "ligand_count": len(ligands),
            "source": "BindingDB best_binding records from IPDB",
            **source_meta,
        },
    }

    for output_path, payload in (
        (args.ipdb_output, ipdb_output),
        (args.ligand_output, ligand_output),
    ):
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")

    print(
        "Wrote "
        f"{len(proteins)} protein names, "
        f"{len(evidence)} evidence records, "
        f"{len(pair_rows)} pairs, and "
        f"{len(ligands)} known-ligand records"
    )
    print(f"IPDB output: {args.ipdb_output}")
    print(f"Ligand output: {args.ligand_output}")


if __name__ == "__main__":
    main()
