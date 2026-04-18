# isovar

**Functions for moving from fine-mapped sQTL variants to isoform-level functional consequences.**

`isovar` is a set of R functions (currently unpackaged, kept in `R/` and loaded via `devtools::load_all()` or `source()`-ing) that implement the variant-annotation and prioritization steps of a larger splicing-genetics workflow aimed at explaining COPD GWAS signals.

## Workflow

```
fine-mapped sQTL variants
        │
        ▼
┌───────────────────────┐
│ (1) prioritize        │  ← splaire / splaireVar deltas, gnomAD AF, rsID
│     variants          │    (isovar)
└───────────┬───────────┘
            ▼
┌───────────────────────┐
│ (2) map variant →     │  ← involved splice site(s) + long-read evidence
│     isoforms          │    (isoscope, extended here first)
└───────────┬───────────┘
            ▼
┌───────────────────────┐
│ (3) functional        │  ← event classification, PTC / NMD / protein
│     consequences      │    (Isopair)
└───────────────────────┘
```

## Inputs

1. Fine-mapped sQTL credible sets (per-variant PIP table).
2. Sequence-based variant-effect splice predictions — `splaire` / `splaireVar` scores.
3. Long-read RNA-seq in the relevant tissue / cell type.
4. GWAS summary statistics (COPD).

## Current scope (what's implemented here now)

| File | Function | Purpose |
|---|---|---|
| `R/splaire_rank.R` | `rankSplaireVariants()` | Rank credible-set variants by predicted splice effect across both model families and three splice-site heads |
| `R/gnomad_annotate.R` | `annotateGnomad()` | Attach rsID + gnomAD v4.1 population allele frequencies via remote tabix over AWS Open Data |
| `R/splice_site_from_splaire.R` | `implicatedSpliceSites()` | Extract structured splice-site objects (chr, pos, head, direction, magnitude, paired-indel-shift) from splaire rows |

Scope will grow toward Step 2 (isoform bridge) once the foundation is exercised against AKR1A1.

## Status

**Unpackaged by design.** Functions are roxygen-documented and tested with `testthat`, but kept as a source-able collection until the API settles. Promotion to a package is a mechanical step when the shape stabilizes.

## Dependencies

See `deps.R`. Currently: `Rsamtools`, `VariantAnnotation`, `GenomicRanges`, `httr2`, `cli`, `cachem`, `memoise`, `readr`, `dplyr`, `tidyr`, `stringr` (all CRAN / Bioconductor).

## Usage

```r
source("deps.R")
for (f in list.files("R", full.names = TRUE)) source(f)

# Prioritize AKR1A1 credible-set variants
ranked <- rankSplaireVariants(
  "HAEC185_sQTL_credible_sets_splaire_scores.tsv.gz",
  gene = "AKR1A1"
)

# Annotate with rsIDs + gnomAD AFs
annotated <- annotateGnomad(ranked$variant_id, dataset = "genomes")

# Derive structured splice-site objects
sites <- implicatedSpliceSites(ranked, magnitude_threshold = 0.1)
```

## Related tools

- [`Isopair`](https://github.com/peter4244/Isopair) — isoform-pair structural event classification (12 event types, PTC/NMD/protein analysis).
- [`isoscope`](https://github.com/peter4244/isoscope) — isoform annotation + long-read extraction + IGV visualization.
- `splaire` / `splaireVar` — sequence-based variant-effect splice prediction (external; outputs consumed here).

## License

MIT.
