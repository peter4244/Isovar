# isovar — function dictionary

One-line purpose per function, the signature, a sentence or two on inputs/outputs, and cross-references. Kept in sync with source; see [`methods.md`](methods.md) for the narrative.

## Variant prioritization & splice-site derivation (`R/splaire_rank.R`, `R/splice_site_from_splaire.R`)

### `rankSplaireVariants(sm_predictions, gene, magnitude_threshold, models, heads)`
Ranks credible-set variants by |delta| across model families (`splaire`, `splaireVar`) × splice-site heads (`don`, `acc`, `ssu`) × direction (`inc`, `dec`).

- `sm_predictions` — path to a splaire TSV or a pre-loaded data.frame (splaire v1 schema).
- `gene` — optional HGNC gene symbol(s) to restrict to.
- `magnitude_threshold` — |delta| bar for the `high_impact` flag (default 0.1).
- **Returns** a tibble with one row per input variant, sorted by `top_abs_delta` desc, carrying variant identity, PIP, and the winning delta's model / head / direction / position / offset.

### `implicatedSpliceSites(sm_predictions, gene, model, heads, magnitude_threshold, motif_shift_threshold, motif_shift_tolerance)`
Long-format splice-site events per variant × head. Classifies each event as `motif_shift_loss` / `motif_shift_gain` / `site_loss` / `site_gain` / `weak`. The motif-shift classification detects indels where the paired gain / loss positions differ by exactly the indel length — evidence the splice motif has been displaced by the indel rather than destroyed or newly created.

- **Returns** a tibble where motif shifts contribute two rows linked by `paired_site_pos`; other events contribute one row.

## Variant annotation (`R/gnomad_annotate.R`)

### `annotateGnomad(variant_ids, dataset, version, endpoint, cache_dir, region_pad_bp)`
Attaches rsID + per-ancestry gnomAD allele frequencies via remote tabix against the AWS Open Data mirror. Region-batched by chromosome; on-disk cache keyed by query hash.

- `variant_ids` — character vector of `chr:pos:ref:alt` IDs.
- `dataset` — `"genomes"` (default, covers intronic variants) or `"exomes"`.
- **Returns** a tibble with one row per input variant (in input order): `rsid`, `filter`, global `af` / `ac` / `an` / `nhomalt`, per-ancestry AFs (`af_afr`, `af_amr`, `af_asj`, `af_eas`, `af_fin`, `af_mid`, `af_nfe`, `af_sas`, `af_remaining`), `grpmax`, `fafmax_faf95_max`.

## GWAS integration (`R/gwas.R`)

### `loadCopdGwas(path, cols)`
Reads a tab-separated GWAS sumstats file with the icgcUkb schema (16 columns). Gzip handled transparently.

### `annotateGwas(variants, gwas, liftover_chain, variants_build, gwas_build)`
Two-pass match of GWAS effect estimates onto a variant table.

- **Pass 1 — rsID** (build-independent, primary).
- **Pass 2 — genomic position** (secondary, for rows that missed rsID). Liftover applied when builds differ; skipped silently if `liftover_chain` is `NULL` and builds differ.
- **Returns** `variants` augmented with `gwas_effectallele`, `gwas_otherallele`, `gwas_eaf`, `gwas_beta`, `gwas_se`, `gwas_p`, `gwas_z`, `gwas_n`, `gwas_imputersq`, `gwas_marker`, `gwas_matched_by` (`"rsid"` / `"position"` / `NA`), `gwas_alignment_status`, `gwas_allele_alignment_ok`.

### `checkAlleleAlignment(df, ref_col, alt_col, effect_col, other_col)`
Classifies ref/alt vs. effect/other orientation as `aligned` / `swapped` / `strand_flip` / `mismatch`. Never silently flips beta — emits a status column and a boolean.

### `getLiftoverChain(from, to, cache_dir)`
Fetches and caches a UCSC liftOver chain file (e.g. `hg38ToHg19`). Returns an absolute path to the uncompressed chain suitable for `rtracklayer::import.chain()` or `annotateGwas(liftover_chain = ...)`.

### `buildAnnotatedCredibleSet(sm_predictions, gwas, gene, gnomad_cache_dir, liftover_chain, variants_build, gwas_build)`
**The Chunk A deliverable.** Composes `rankSplaireVariants` → `annotateGnomad` → `annotateGwas` into a single merged per-variant table. One function, one tibble, one gene.

## Conventions

- **Variant IDs** are `chr:pos:ref:alt` (with or without the `chr` prefix; normalization happens at the boundary).
- **Build defaults** — variants (splaire, gnomAD) GRCh38; GWAS GRCh37. Always overridable.
- **Cache directories** — all remote-data functions accept a `cache_dir` argument; none are required. Default behavior is no cache.
- **Argument name `sm_predictions`** — generic for splicing-model predictions. Currently splaire v1 schema only; see `docs/methods.md` "Design notes".

## Pending (chunks B–G)

Not yet implemented:

- `groupHaplotypes()` (Chunk B)
- `runIsoscopeGene()`, `classifyIsoformsBySiteUsage()` (Chunk C)
- `longreadSources()`, `loadLongreadEvidence()` (Chunk D)
- `parseStructuresMultiGtf()`, `loadTranscriptSequences()`, `hiddenPtcStatus()`, `computeDominantIsoform()` (Chunk E)
- `buildSpliceReport()` (Chunk F)
