# isovar — data schemas

Explicit schemas for every input and output isovar touches. Each schema lists **required** columns (absence is an error) and **optional** columns (preserved when present). A column glossary at the end defines any fields whose meaning isn't obvious.

This document is authoritative — functions validate against these schemas at boundaries.

## Provenance metadata

Every isovar result carries a structured provenance metadata object, attached as an R attribute (`attr(x, "isovar_meta")`) or written as a sidecar `.meta.json` next to any saved file.

```
list(
  isovar_version = "0.0.2",                  # isovar version or git SHA
  generated_at   = "2026-04-18T11:30:00Z",   # ISO 8601 timestamp
  genome_build   = "GRCh38",                 # object's coordinate build
  gtf_version    = "GENCODE_v49",            # annotation version when applicable
  sources        = list(                     # list of upstream inputs
    list(
      id      = "haec185_splaire",
      kind    = "splicing_model_predictions",
      path    = "/.../HAEC185_sQTL_credible_sets_splaire_scores.tsv.gz",
      model   = "splaire / splaireVar (v1)",
      cohort  = "HAEC-185",
      tissue  = "basal airway epithelium",
      genome_build = "GRCh38",
      hash    = NA_character_            # optional sha256
    ),
    list(
      id      = "gnomad_v4.1_genomes",
      kind    = "variant_annotation",
      endpoint = "https://gnomad-public-us-east-1.s3.amazonaws.com",
      version = "4.1",
      dataset = "genomes",
      genome_build = "GRCh38"
    ),
    list(
      id      = "icgcUkb_copd_gwas",
      kind    = "gwas_sumstats",
      path    = "/.../icgcUkb-20251213.dn8.gz",
      phenotype = "COPD",
      genome_build = "GRCh37"
    )
  )
)
```

**Merge rule:** when multiple sources are combined, the output's `sources` list is the concatenation of the inputs' `sources` — never collapsed. `genome_build` / `gtf_version` are carried from the primary (left) object and flagged if inputs disagree.

**Accessors:**
- `getIsovarMeta(x)` — returns the metadata list.
- `setIsovarMeta(x, ...)` — creates or updates the attribute.
- `mergeIsovarMeta(...)` — merges metadata from multiple objects, preserving all source records.
- `writeMeta(x, path)` — writes a sidecar JSON.
- `readMeta(path)` — reads a sidecar JSON.

## Input schemas

### Splicing-model predictions (`sm_predictions`)

**Current accepted format: splaire v1 TSV** (the 80-column table produced for HAEC-185). Validation is against column names; values are not type-checked at the boundary beyond what `readr` infers.

Required:
- `gene` (HGNC symbol)
- `phenotype` (intron-cluster id, `chr:start:end:clu_N_+/-`)
- `ensembl_id`
- `variant_id` (`chr{N}:{pos}:{ref}:{alt}`)
- `credible_set_number` (integer)
- `posterior_inclusion_probability` (numeric in `[0, 1]`)
- `chr`, `strand`
- For each `{splaire, splaireVar} × {don, acc, ssu} × {max_inc, max_dec}`:
  - the value (`splaire_don_max_inc`, …)
  - `_off` (offset from variant in bp)
  - `_pos` (genomic position of the max)
  - — 36 columns total
- For each `{splaire, splaireVar} × {don, acc, ssu} × {istart, iend} × {ref, alt, delta}`: 36 columns. `NA` when the junction falls outside the ±5 kb scored window.

Planned future: a normalized generic schema that any splicing model can coerce to, with per-source metadata carried separately.

### Variant annotation (gnomAD output, `annotateGnomad` return)

**Default fields** (returned unless `fields = "all"`):
- `variant_id`, `chr`, `pos`, `ref`, `alt`
- `rsid`, `filter`
- `af`, `af_nfe`, `af_eas`, `grpmax`

**Extended fields** (only with `fields = "all"`):
- `ac`, `an`, `nhomalt`
- `af_afr`, `af_amr`, `af_asj`, `af_fin`, `af_mid`, `af_sas`, `af_remaining`
- `fafmax_faf95_max`, `fafmax_faf95_max_gen_anc`

### GWAS sumstats (`loadCopdGwas` return, icgcUkb schema)

Required columns in the source file:

- `marker` (`chr:pos:ref:alt`, no `chr` prefix)
- `chr`, `pos`, `rsid`
- `effectallele`, `otherallele`
- `eaf`, `beta`, `se`, `p`
- `z`, `imputersq`, `n`
- `pheno`, `mode`, `tag`

## Output schemas

### Annotated credible-set (primary Chunk A output)

**Shape:** a nested tibble — one row per credible set, with a list-column `variants` holding the per-variant tibble. Default shape returned by `buildAnnotatedCredibleSet()`. Helper `as_flat(x)` unnests into the legacy flat form.

**Credible-set-level row:**

- `credible_set_id` — synthesized as `paste(gene, phenotype, sep = "__")`. Unique within a single build/run.
- `gene` (HGNC symbol)
- `phenotype` (intron cluster id)
- `ensembl_id`
- `strand`
- `credible_set_number` (integer)
- `n_variants` (integer)
- `variants` (list-column; see below)

**Per-variant tibble inside `variants`:**

Identity / splaire:
- `variant_id`, `chr`, `pos`, `ref`, `alt`, `rsid`
- `pip` (renamed from `posterior_inclusion_probability`)
- `top_model`, `top_head`, `top_direction`, `top_delta`, `top_abs_delta`, `top_pos`, `top_off`, `high_impact`

gnomAD (default columns):
- `filter`, `af`, `af_nfe`, `af_eas`, `grpmax`

GWAS (default columns):
- `gwas_effectallele`, `gwas_otherallele`, `gwas_eaf`
- `gwas_beta`, `gwas_se`, `gwas_p`
- `gwas_matched_by` (`"rsid"` / `"position"` / `NA`)
- `gwas_alignment_status` (see glossary)
- `gwas_allele_alignment_ok` (logical)

### Flat form (output of `as_flat()`)

The nested output unnested into a single flat tibble — one row per variant — with credible-set-level columns duplicated across the CS's variant rows. Column set = credible-set columns + per-variant columns above. Retained for ad-hoc viewing, not as the canonical shape.

### Disk format (output of `write_tsv_pair()`)

Two TSVs in a target directory:
- `credible_sets.tsv` — the credible-set-level rows, minus the `variants` list-column.
- `variants.tsv` — flattened per-variant rows with a leading `credible_set_id` column joining back to `credible_sets.tsv`.

Both TSVs get a `.meta.json` sidecar with the provenance metadata.

## Column glossary

- **`pip`** (or `posterior_inclusion_probability`) — fine-mapping posterior inclusion probability for the variant within its credible set. In `[0, 1]`; sums over a credible set at or near 1 by construction.
- **`top_abs_delta`** — largest absolute splice-prediction delta (alt − ref) across the configured model families and splice-site heads in a ±5 kb scanning window around the variant. Often abbreviated `|Δ|` in narrative. Bounded on `[0, 1]`. **Values near the ceiling (~0.85+) are in the model's saturation regime — small numerical differences between such variants should not be treated as meaningful causal rankings.** Tie-break with phenotype alignment (which splice site does the variant affect, relative to the measured phenotype?), long-read observations, or conditional analysis, not with `top_abs_delta` precision.
- **`top_model`** — which model family produced the winning delta: `splaire` (reference-trained) or `splaireVar` (variant-trained).
- **`top_head`** — which splice-site head: `don` (donor / 5′ss), `acc` (acceptor / 3′ss), `ssu` (general splice-site usage).
- **`top_direction`** — sign of the winning delta: `inc` (creation / strengthening) or `dec` (loss / weakening).
- **`high_impact`** — boolean, `top_abs_delta >= magnitude_threshold` (default 0.1).
- **`filter`** — gnomAD site-level filter (`PASS` or a comma-separated list of failure codes).
- **`af`** — gnomAD global allele frequency.
- **`af_nfe`, `af_eas`, …** — gnomAD per-ancestry allele frequencies.
- **`grpmax`** — **gnomAD "group max"**: the highest allele frequency across all reported ancestries, labeled with which ancestry had it (e.g. `grpmax = "eas"` when `af_eas` is the max). Useful as a rarity screen — `grpmax > 0.05` means the variant is common in at least one ancestry; `grpmax < 0.001` means rare everywhere.
- **`gwas_matched_by`** — `"rsid"` if the GWAS row was joined by rsID, `"position"` if by genomic coordinate (with optional liftover), `NA` if no match.
- **`gwas_alignment_status`** — how the variant's `ref/alt` relate to the GWAS's `effectallele/otherallele`:
  - `aligned` — `alt == effectallele` and `ref == otherallele`. β refers to the alt allele, as expected.
  - `swapped` — `ref == effectallele` and `alt == otherallele`. β refers to the *reference* allele; flip the sign if using.
  - `strand_flip` — alleles match only after reverse-complementing both sides. Ambiguous and not claimed for palindromic (A/T or C/G) SNPs.
  - `mismatch` — no orientation works. Do not use the β.
  - `NA` — no GWAS record matched.
- **`gwas_allele_alignment_ok`** — `TRUE` only when `gwas_alignment_status == "aligned"`.
