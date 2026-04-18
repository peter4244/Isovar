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

### `buildAnnotatedCredibleSet(sm_predictions, gwas, gene, gnomad_cache_dir, liftover_chain, variants_build, gwas_build, fields)`
**The Chunk A deliverable.** Composes `rankSplaireVariants` → `annotateGnomad` → `annotateGwas` and returns a **nested credible-set tibble** — one row per credible set, with a list-column `variants` holding the per-variant tibble. `fields = "default"` (the default) prunes gnomAD and GWAS columns to the canonical set documented in `docs/schemas.md`; `fields = "all"` carries every upstream column forward. Provenance metadata from all three sources is attached via [setIsovarMeta()].

## Output shape helpers (`R/gwas.R`)

### `nestCredibleSets(flat)`
Re-shapes a flat per-variant tibble (from the three composed annotation functions) into the canonical nested credible-set form.

### `as_flat(x)`
Inverse of `nestCredibleSets()`. Unnests the `variants` list-column back into a single row per variant, duplicating credible-set-level columns.

### `write_tsv_pair(x, dir)`
Writes a nested credible-set tibble as two TSVs — `credible_sets.tsv` (one row per CS) + `variants.tsv` (one row per variant with `credible_set_id` key) — plus `.meta.json` sidecars next to each.

## Provenance metadata (`R/metadata.R`)

Every isovar-produced object carries a structured `isovar_meta` attribute. See `docs/schemas.md` for the full schema.

### `setIsovarMeta(x, meta, ...)`
Attach / overwrite metadata. Accepts either a full metadata list (`setIsovarMeta(x, merged_list)`) or named fields (`setIsovarMeta(x, genome_build = "GRCh38", sources = list(...))`).

### `getIsovarMeta(x, require)`
Retrieve the attached metadata. Errors by default when absent; pass `require = FALSE` for a `NULL` return instead.

### `mergeIsovarMeta(...)`
Concatenate metadata from multiple objects, dedup'd by source `id`. Flags `mixed_builds` / `mixed_gtf` when inputs disagree (does not silently resolve).

### `writeMeta(x, data_path)` / `readMeta(data_path)`
Sidecar JSON round-trip. Writes `<data_path>.meta.json`; reads it back on load.

## Conventions

- **Variant IDs** are `chr:pos:ref:alt` (with or without the `chr` prefix; normalization happens at the boundary).
- **Build defaults** — variants (splaire, gnomAD) GRCh38; GWAS GRCh37. Always overridable.
- **Cache directories** — all remote-data functions accept a `cache_dir` argument; none are required. Default behavior is no cache.
- **Argument name `sm_predictions`** — generic for splicing-model predictions. Currently splaire v1 schema only; see `docs/methods.md` "Design notes".

## Haplotype grouping (`R/haplotype.R`)

### `groupHaplotypes(variants, population, r2_threshold, token, fixture)`
Groups variants into haplotypes using LDlink r² (via `LDlinkR::LDmatrix`) with hierarchical clustering at `r2_threshold` (default 0.8). Accepts flat per-variant tibbles or nested credible-set tibbles. For variants that aren't in the 1000 Genomes panel (most indels), falls back to **proximity assignment** — the orphan inherits the haplotype of its nearest-position neighbour. The `haplotype_assigned_by` column surfaces `"ld_cluster"` vs. `"proximity"` so downstream code can treat proximity calls with appropriate skepticism.

Adds four columns:
- `haplotype_id` (`"hap_1"`, `"hap_2"`, …)
- `causal_candidate` (logical; TRUE for the variant with the largest `top_abs_delta` within each haplotype)
- `ld_r2_to_causal` (numeric; r² to this haplotype's causal candidate)
- `haplotype_assigned_by` (`"ld_cluster"` or `"proximity"`)

Requires `Sys.getenv("LDLINK_TOKEN")` for the live path; pass `fixture = <RDS path>` for an offline / CI path. See `scripts/build_haplotype_fixture.R` to generate a fixture.

## Pending (chunks C–G)

Not yet implemented:

- `runIsoscopeGene()`, `classifyIsoformsBySiteUsage()` (Chunk C)
- `longreadSources()`, `loadLongreadEvidence()` (Chunk D)
- `parseStructuresMultiGtf()`, `loadTranscriptSequences()`, `hiddenPtcStatus()`, `computeDominantIsoform()` (Chunk E)
- `buildSpliceReport()` (Chunk F)
