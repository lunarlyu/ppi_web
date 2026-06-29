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

Regenerate the STRING cache with:

```bash
python3 scripts/fetch_string_cache.py TP53 EGFR BRCA1 MYC AKT1 MTOR
```

Process DepMap Chronos/CERES gene dependency scores with:

```bash
python3 scripts/process_depmap_essential.py CRISPRGeneEffect.csv \
  --score-threshold -0.5 \
  --min-dependency-fraction 0.5
```

DepMap gene effect values are dependency scores, so a lower score means a
stronger dependency. The default rule marks a gene as essential when it scores
`<= -0.5` in at least half of the profiled cell lines with at least 50 valid
scores.

The default STRING confidence threshold is `0.99`, which maps to
`required_score=990` in the STRING API. That is intentionally stringent so
highly connected proteins such as `TP53` open with a readable star map; the
slider can be relaxed toward STRING's broader high-confidence `0.70` cutoff.
