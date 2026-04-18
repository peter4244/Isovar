# isovar — methods

Evolving methods document for the isovar R codebase. Updated per chunk. For the function-level reference see [`functions.md`](functions.md).

## Goal

Produce a structured per-(gene, haplotype) report that moves from **fine-mapped sQTL variants** to **isoform-level functional consequence** — in a form that supports colocalization and interpretation against COPD GWAS. isovar is a thin orchestration layer that wires together three upstream tools:

- **splaire / splaireVar** — sequence-based variant-effect splice-site predictions (external). isovar treats these as "splicing-model predictions" (`sm_predictions` in function signatures) because we expect the input to grow beyond splaire: multi-model, multi-tissue predictions are a planned extension.
- **Isopair** — isoform-pair structural event classification + PTC/NMD/protein consequences. R package, library dependency.
- **isoscope** — GENCODE + SQANTI isoform annotation + IGV visualization. Subprocess call to the canonical checkout.

Code-reuse rule: never duplicate logic from Isopair or isoscope. Fix bugs at their source; isovar should stay a thin bridge.

## Stages

### Stage 1 — variant prioritization and merge (Chunk A)

One composed call — `buildAnnotatedCredibleSet()` — produces a per-variant tibble combining:

1. **Fine-mapped sQTL credible-set rows** — carried through from the input splicing-model table (`phenotype`, `credible_set_number`, `posterior_inclusion_probability`, `ensembl_id`, per-variant PIP).
2. **Splicing-model deltas** — for each variant, the winning `|delta|` across two model families (`splaire`, `splaireVar`) and three splice-site heads (donor `don`, acceptor `acc`, splice-site usage `ssu`), with direction (`inc` / `dec`), position, and offset from the variant. Implicated splice sites are separately derivable via `implicatedSpliceSites()` — including detection of the "motif shifted by indel length" pattern (paired `max_inc` / `max_dec` at an offset exactly equal to an indel's length).
3. **gnomAD annotation** — rsID + per-ancestry allele frequencies via remote tabix against the AWS Open Data gnomAD v4.1 mirror (free egress). Region-batched with an on-disk cache.
4. **GWAS effect estimates** — two-pass match:
   - **Pass 1: rsID** (build-independent, primary).
   - **Pass 2: genomic position** (secondary, requires a liftover chain when variants and GWAS are on different builds). A `gwas_matched_by` column surfaces which pass resolved each row.

Allele alignment between variant ref/alt and GWAS effect/other is classified as `aligned` / `swapped` / `strand_flip` / `mismatch` — beta is never silently flipped.

### Stage 2 — isoform bridge (Chunks C–E, planned)

For each implicated splice site, annotate the gene's isoforms from multiple long-read sources (HAEC-185 basal + NMD lungcells, with cell-type / treatment column masks) and classify which isoforms use canonical vs. shifted splice sites. Identify the dataset-specific dominant isoform and pair it against each comparator. Feed pairs to Isopair for event classification, PTC/hidden-PTC (via `traceReferenceAtg`), NMD attribution, and protein consequence. isoscope is invoked as a subprocess per source.

### Stage 3 — haplotype grouping (Chunk B, planned)

Variants are grouped into haplotypes via **LDlinkR** (official CRAN wrapper around the LDlink REST API), using actual r² / D′ rather than AF-cluster heuristics. The splaire-prioritized causal candidate is flagged within each haplotype. Tag SNPs that hitchhike with an indel causal variant get correctly attributed to their haplotype — something fine-mapping alone cannot do.

### Stage 4 — per (gene, haplotype) report (Chunk F, planned)

A tibble-of-lists shape per (gene, haplotype) with fields for variants, implicated sites (with per-LR-source + per-condition observation flags), isoform pairs, events, PTC status, hidden-PTC status, protein effects, and GWAS association. Rendering (Rmd HTML, IGV sessions) is a separate concern handled later.

## Design decisions (locked in)

- **Genome build / annotation**: GRCh38 for splaire / gnomAD / HAEC-185; GRCh37 for the COPD GWAS. Defaults follow that but are always explicit parameters; per-project version metadata is tracked rather than assumed.
- **GENCODE v49** is the default annotation version.
- **rsID-first, position-second** for GWAS matching. Position matching runs routinely on every call, with optional liftover for cross-build joins. For the icgcUkb GWAS specifically, position matching adds zero matches because the GWAS systematically excludes indels — this is surfaced by the `gwas_matched_by` column (`rsid` / `position` / `NA`) and by the absence of indel records entirely.
- **Haplotype interpretation is first-class.** GWAS significance on tag SNPs alone cannot identify indel causal variants that aren't in the sumstats; LDlink-based grouping is what ties them together.
- **Long-read evidence is qualitative** at our current sample size (9 HAEC-185 basal donors + the NMD cohort). Quantitative allele effects come from short-read leafcutter sQTL PSI; long-read data confirms isoform structure and supports Isopair event classification.
- **NMD analysis is systematic** — `computePtcStatus` + `traceReferenceAtg` are applied to every coding or non-coding-but-CDS-tractable isoform pair, regardless of event type. Novel intron-retention isoforms (common in SQANTI output) are PTC-invisible under `computePtcStatus` alone; `traceReferenceAtg` recovers the hidden cases by tracing the reference ATG through the comparator reading frame.

## Design notes: multi-source / multi-tissue splicing predictions

The `sm_predictions` argument is named generically in anticipation of:

- Multiple fine-mapped sQTL sources (e.g. HAEC-185 basal, ALI, GTEx-Lung, GTEx-Airway).
- Multiple splicing-effect models (splaire, splaireVar, SpliceAI, Pangolin, tissue-retrained variants).
- Per-source metadata (tissue, cohort, model version, assembly, training corpus) that rides with the predictions.

The current schema requirement is splaire v1 (the 80-column TSV produced for HAEC-185). The eventual shape is a list of source objects, each carrying `{id, tissue, model, schema, data}`, with downstream code summarizing "most-significant effect per variant across sources" and flagging discordance between sources.

Until that refactor lands, the documented guarantee is: the function accepts splaire v1 only; it errors explicitly on unrecognized schemas rather than guessing.

## Provenance metadata

Every isovar-produced object carries an `isovar_meta` attribute recording: isovar version, generation timestamp, genome build, GTF/annotation version (when applicable), and a list of upstream `sources` (each with `id`, `kind`, path/endpoint, model/cohort/tissue/build). Metadata propagates through every merge — when multiple inputs combine, their sources concatenate rather than collapse, and build/GTF disagreements are flagged as `mixed_builds` / `mixed_gtf` fields rather than silently resolved. On write, metadata is emitted as a sidecar `.meta.json` next to every saved TSV/RDS.

See `docs/schemas.md` ("Provenance metadata") for the full schema and `R/metadata.R` for the accessors (`setIsovarMeta`, `getIsovarMeta`, `mergeIsovarMeta`, `writeMeta`, `readMeta`).

## Change log

- 2026-04-18 — Chunk A landed: `buildAnnotatedCredibleSet()`, `annotateGwas()` (rsID + position match), `loadCopdGwas()`, `checkAlleleAlignment()`, `getLiftoverChain()`. Generic `sm_predictions` naming adopted throughout.
- 2026-04-18 — Chunk A.3 landed: explicit input/output schemas (`docs/schemas.md`), nested credible-set output shape (one row per CS with list-column `variants`), helpers `nestCredibleSets()` / `as_flat()` / `write_tsv_pair()`, default column pruning for gnomAD + GWAS (opt in via `fields = "all"`), and provenance metadata propagation (`R/metadata.R`, sidecar JSON on write).
