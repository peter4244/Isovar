# isovar

## 0.0.2 (unreleased)

- `loadCopdGwas()` — reads icgcUkb COPD summary statistics.
- `annotateGwas()` — two-pass match: rsID (build-independent, primary) then genomic position with optional liftover (secondary). `gwas_matched_by` column surfaces which pass resolved each row.
- `checkAlleleAlignment()` — classifies effect-allele orientation (`aligned` / `swapped` / `strand_flip` / `mismatch`) without silently flipping beta.
- `getLiftoverChain()` — fetches + caches UCSC chain files (default: hg38ToHg19).
- `buildAnnotatedCredibleSet()` — composed Step-1 pipeline returning a **nested credible-set tibble** (one row per CS with list-column `variants`).
- `nestCredibleSets()`, `as_flat()`, `write_tsv_pair()` — output-shape helpers.
- `setIsovarMeta()`, `getIsovarMeta()`, `mergeIsovarMeta()`, `writeMeta()`, `readMeta()` — provenance-metadata propagation; every isovar result carries genome build + GTF version + source list + timestamp; sidecar `.meta.json` on write.
- `sm_predictions` is the canonical generic argument name for splicing-model variant-effect predictions (splaire v1 schema for now; multi-model / multi-tissue planned).
- Default column pruning: gnomAD returns `rsid, filter, af, af_nfe, af_eas, grpmax` by default (`fields = "all"` for the extended set); GWAS drops `gwas_z`/`gwas_n`/`gwas_imputersq`/`gwas_marker` from the default merged table.
- Explicit schemas + column glossary documented in `docs/schemas.md`. Methods narrative in `docs/methods.md`; function dictionary in `docs/functions.md`.
- `groupHaplotypes()` — haplotype grouping via `LDlinkR::LDmatrix` + hierarchical clustering at `r2_threshold`. Proximity fallback for variants absent from 1KG (most indels). Adds `haplotype_id`, `causal_candidate`, `ld_r2_to_causal`, `haplotype_assigned_by` columns. Live path gated by `LDLINK_TOKEN`; CI uses an RDS fixture.
- `scripts/build_haplotype_fixture.R` — one-time helper to pre-compute the AKR1A1 r² matrix for offline tests.
- `loadIsovarSecrets()` — user-level `KEY=VALUE` secrets at `~/.config/isovar/secrets.env`.
- `runIsoscopeGene()` — thin wrapper around isoscope's `gene_isoform_annotation.R` with `--config <path>` (upstream isoscope change); each long-read source drives isoscope through its own `inst/isoscope_configs/*.R` file. `ISOVAR_SOURCE_META` in each config propagates into the output's `isovar_meta`.
- `classifyIsoformsBySiteUsage()` — per-isoform logical columns for whether a junction is anchored at each requested genomic site.
- `inst/isoscope_configs/haec185_basal.R` and `inst/isoscope_configs/nmd_lungcells.R` — canonical source configs for the two long-read datasets in scope.
- Optional **isocall (pre-SQANTI) mode** via `inst/isoscope_configs/{haec185_basal,nmd_lungcells}_isocall.R`. Requires the isoscope upstream change adding `SOURCE_KIND = "isocall"` (isoscope@157f121). Useful for testing whether SQANTI's junction correction has altered the long-read evidence relative to raw calls.
- `longreadSources()`, `parseSampleColumns()`, `loadLongreadEvidence()` — multi-source long-read evidence with per-condition stratification. Per-source count matrices are parsed into donor + celltype + treatment metadata via the source's `sample_schema` (`donor_only` for HAEC-185, `sample_celltype_donor_treatment` for NMD). NMD sources honor `celltypes` / `treatments` filters at load time.
- `loadSqtlResults()` / `annotateSqtl()` — short-read Leafcutter sQTL ingestion (tensorqtl `cis_qtl_pairs` parquet format) with arrow-backed filter pushdown. `inst/sqtl_configs/haec185_basal_leafcutter.R` is the canonical config.

## 0.0.1 (unreleased)

Initial scaffold.

- `rankSplaireVariants()` — rank credible-set variants by splaire / splaireVar delta across donor / acceptor / splice-site-usage heads.
- `annotateGnomad()` — attach rsID and gnomAD v4.1 population allele frequencies via remote tabix over AWS Open Data (region-batched).
- `implicatedSpliceSites()` — derive structured splice-site objects from splaire rows, detecting the "motif shifted by indel length" pattern.
- First runs prototyped against AKR1A1 credible-set variants from HAEC-185.
