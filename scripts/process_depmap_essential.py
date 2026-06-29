#!/usr/bin/env python3
"""Convert DepMap Chronos/CERES gene dependency scores into browser JSON."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import DefaultDict, TextIO


GENE_WITH_ENTREZ_RE = re.compile(r"^(.+?)\s+\(\d+\)$")


def clean_symbol(column_name: str) -> str:
    """Convert DepMap headers such as 'TP53 (7157)' to 'TP53'."""
    name = column_name.strip()
    match = GENE_WITH_ENTREZ_RE.match(name)
    if match:
        name = match.group(1)
    return name.strip().upper()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def parse_score(value: str) -> float | None:
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None

    return score if math.isfinite(score) else None


def infer_essentials(
    input_path: Path,
    score_threshold: float,
    min_dependency_fraction: float,
    min_cell_lines: int,
) -> tuple[list[str], dict[str, dict[str, float | int]]]:
    stats: DefaultDict[str, dict[str, float | int]] = defaultdict(
        lambda: {"dependent": 0, "score_sum": 0.0, "valid": 0}
    )

    with open_text(input_path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        symbols = [clean_symbol(column) for column in header[1:]]

        for row in reader:
            for offset, symbol in enumerate(symbols, start=1):
                if offset >= len(row) or not symbol:
                    continue

                score = parse_score(row[offset])
                if score is None:
                    continue

                gene_stats = stats[symbol]
                gene_stats["valid"] += 1
                gene_stats["score_sum"] += score
                if score <= score_threshold:
                    gene_stats["dependent"] += 1

    summary = {}
    essential_symbols = []

    for symbol, gene_stats in stats.items():
        valid = int(gene_stats["valid"])
        if valid < min_cell_lines:
            continue

        dependent = int(gene_stats["dependent"])
        fraction = dependent / valid
        summary[symbol] = {
            "dependent_cell_lines": dependent,
            "dependency_fraction": round(fraction, 6),
            "mean_score": round(float(gene_stats["score_sum"]) / valid, 6),
            "profiled_cell_lines": valid,
        }

        if fraction >= min_dependency_fraction:
            essential_symbols.append(symbol)

    return sorted(essential_symbols), dict(sorted(summary.items()))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Infer an essential protein list from a DepMap Chronos/CERES gene "
            "dependency matrix. Scores <= --score-threshold are counted as "
            "dependencies per cell line."
        )
    )
    parser.add_argument("input", type=Path, help="DepMap gene effect CSV or CSV.GZ")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/essential_proteins.json"),
        help="Processed JSON consumed by the browser app",
    )
    parser.add_argument(
        "--score-threshold",
        type=float,
        default=-0.5,
        help="A gene counts as dependent in one cell line at or below this score",
    )
    parser.add_argument(
        "--min-dependency-fraction",
        type=float,
        default=0.5,
        help="Minimum fraction of profiled cell lines where the gene is dependent",
    )
    parser.add_argument(
        "--min-cell-lines",
        type=int,
        default=50,
        help="Minimum non-empty scores required for a gene to be considered",
    )
    args = parser.parse_args()

    if not 0 <= args.min_dependency_fraction <= 1:
        raise ValueError("--min-dependency-fraction must be between 0 and 1")

    proteins, summary = infer_essentials(
        args.input,
        args.score_threshold,
        args.min_dependency_fraction,
        args.min_cell_lines,
    )
    mean_scores = [
        float(gene_summary["mean_score"]) for gene_summary in summary.values()
    ]
    payload = {
        "meta": {
            "description": (
                "Gene-level summaries inferred from DepMap Chronos/CERES "
                "dependency scores, plus an essential-protein list."
            ),
            "essential_gene_count": len(proteins),
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "mean_score_range": {
                "max": max(mean_scores) if mean_scores else None,
                "min": min(mean_scores) if mean_scores else None,
            },
            "min_cell_lines": args.min_cell_lines,
            "min_dependency_fraction": args.min_dependency_fraction,
            "profiled_gene_count": len(summary),
            "score_threshold": args.score_threshold,
            "source": "DepMap Chronos/CERES gene dependency scores",
            "source_file": args.input.name,
            "symbol_case": "uppercase",
        },
        "proteins": proteins,
        "summary": summary,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        f"Wrote {len(proteins)} essential proteins and "
        f"{len(summary)} scored genes to {args.output}"
    )


if __name__ == "__main__":
    main()
