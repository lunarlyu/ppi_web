# STRING PPI Star Map

An interactive, static STRING protein-protein interaction explorer.

## Run locally

```bash
python3 -m http.server 8000
```

Then open `http://127.0.0.1:8000`.

## Data

- `data/string_interactions.json` is a processed cache for search suggestions and
  offline fallback. The app still queries STRING directly for the current
  protein so searches are not limited to the bundled cache.
- `data/essential_proteins.json` is the processed local DepMap dependency file.
  It stores the inferred essential-protein list in `proteins` and mean Chronos /
  CERES dependency statistics for every profiled gene in `summary`.
- `data/ipdb_processed.json` is the compact IPDB extract for RIPTAC pairs and
  protein-card evidence. It is generated from the upstream `ipdb-snapshot.json`
  in `~/Downloads/Induced Proximity Database Project copy/public/`.
- `data/known_ligands.json` is the BindingDB best-binding index generated from
  the same upstream IPDB snapshot.

Regenerate the STRING cache with:

```bash
python3 scripts/fetch_string_cache.py TP53 EGFR BRCA1 MYC AKT1 MTOR
```

Process DepMap Chronos/CERES gene dependency scores with:

```bash
python3 scripts/process_depmap_essential.py CRISPRGeneEffect.csv \
  --score-threshold -1.0 \
  --min-dependency-fraction 0.5
```

DepMap gene effect values are dependency scores, so a lower score means a
stronger dependency. The default rule marks a gene as essential when it scores
`<= -1.0` in at least half of the profiled cell lines with at least 50 valid
scores.

Process IPDB evidence and BindingDB best bindings with:

```bash
python3 scripts/process_ipdb_snapshot.py
```

The default STRING confidence threshold is `0.99`, which maps to
`required_score=990` in the STRING API. That is intentionally stringent so
highly connected proteins such as `TP53` open with a readable star map; the
slider can be relaxed toward STRING's broader high-confidence `0.70` cutoff.
