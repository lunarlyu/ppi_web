#!/usr/bin/env python3
"""Fetch a compact Homo sapiens STRING cache for the static PPI star map."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


DEFAULT_PROTEINS = ["TP53", "EGFR", "BRCA1", "MYC", "AKT1", "MTOR"]
EVIDENCE_FIELDS = ["nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"]
API_URL = "https://string-db.org/api/tsv/interaction_partners"


def fetch_interactions(protein: str, species: int, required_score: int, limit: int) -> str:
    params = urllib.parse.urlencode(
        {
            "identifier": protein,
            "species": species,
            "required_score": required_score,
            "limit": limit,
            "caller_identity": "ppi-star-map-cache",
        }
    )
    request = urllib.request.Request(
        f"{API_URL}?{params}",
        headers={"User-Agent": "ppi-star-map-cache/1.0"},
    )

    with urllib.request.urlopen(request, timeout=20) as response:
        return response.read().decode("utf-8")


def parse_tsv(tsv: str) -> list[dict]:
    rows = []
    for row in csv.DictReader(tsv.splitlines(), delimiter="\t"):
        rows.append(
            {
                "stringId_A": row["stringId_A"],
                "stringId_B": row["stringId_B"],
                "preferredName_A": row["preferredName_A"],
                "preferredName_B": row["preferredName_B"],
                "score": float(row["score"]),
                "evidence": {
                    field: float(row.get(field) or 0)
                    for field in EVIDENCE_FIELDS
                },
            }
        )
    return rows


def dedupe(rows: list[dict]) -> list[dict]:
    seen: dict[tuple[str, str], dict] = {}

    for row in rows:
        key = tuple(sorted((row["stringId_A"], row["stringId_B"])))
        if key not in seen or row["score"] > seen[key]["score"]:
            seen[key] = row

    return sorted(
        seen.values(),
        key=lambda row: (-row["score"], row["preferredName_A"], row["preferredName_B"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("proteins", nargs="*", default=DEFAULT_PROTEINS)
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument("--required-score", type=int, default=700)
    parser.add_argument("--limit", type=int, default=80)
    parser.add_argument("--output", type=Path, default=Path("data/string_interactions.json"))
    args = parser.parse_args()

    interactions = []
    for index, protein in enumerate(args.proteins):
        if index:
            time.sleep(0.25)
        interactions.extend(
            parse_tsv(
                fetch_interactions(
                    protein.upper(),
                    args.species,
                    args.required_score,
                    args.limit,
                )
            )
        )

    payload = {
        "meta": {
            "description": "Processed STRING interaction cache for the static PPI star map.",
            "species": args.species,
            "required_score": args.required_score,
            "limit_per_seed": args.limit,
            "seed_proteins": [protein.upper() for protein in args.proteins],
            "source": "STRING interaction_partners API",
        },
        "interactions": dedupe(interactions),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {len(payload['interactions'])} interactions to {args.output}")


if __name__ == "__main__":
    main()
