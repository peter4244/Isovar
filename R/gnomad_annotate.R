#' Annotate variants with rsID and gnomAD population allele frequencies
#'
#' Given a character vector of `chr:pos:ref:alt` variant IDs (GRCh38),
#' fetch the matching records from gnomAD v4.1 sites VCFs via remote
#' tabix over HTTPS (default endpoint: the AWS Open Data mirror, which
#' is free to read — see README). Returns one row per input variant,
#' preserving input order; variants absent from gnomAD get `NA` values
#' in the gnomAD-derived columns.
#'
#' Queries are region-batched per chromosome: variants grouped by
#' `chr`, one `scanTabix()` call pulls the enclosing range, and the
#' returned records are filtered to exact `chr:pos:ref:alt` matches.
#' This is ~5–10× faster than per-variant queries for clustered inputs
#' (credible sets, exon regions).
#'
#' Multi-allelic sites are handled: matches are on full
#' `chr:pos:ref:alt`, not just position. If a site has multiple ALTs in
#' gnomAD, only the one matching the input's ALT is returned.
#'
#' @param variant_ids Character vector of IDs formatted `chr:pos:ref:alt`
#'   (e.g. `"chr1:45549863:A:ATATCTG"`). A leading `chr` prefix is
#'   accepted or omitted; internally normalized to match gnomAD v4 style.
#' @param dataset `"genomes"` (default) or `"exomes"`. Use `"genomes"`
#'   for intronic / regulatory variants; `"exomes"` only covers coding
#'   + flanking regions.
#' @param version gnomAD release version. Default `"4.1"`. Other values
#'   assume the same URL template.
#' @param endpoint Base URL up to but not including the dataset path.
#'   Default: AWS Open Data (`https://gnomad-public-us-east-1.s3.amazonaws.com`).
#'   Configurable to point at Azure/GCP mirrors or a local copy.
#' @param cache_dir Optional directory for on-disk caching of parsed
#'   region pulls. If supplied, repeated queries for overlapping
#'   regions skip the network. Default `NULL` (no caching).
#' @param region_pad_bp Integer; widen each region's query span by this
#'   many base pairs on each side. Default 0. Raise if you expect
#'   nearby variants to be queried in follow-up calls (helps cache hits).
#' @return A tibble with one row per input `variant_ids` entry, in
#'   input order, with columns:
#'   - `variant_id`, `chr`, `pos`, `ref`, `alt` — input echoes
#'   - `rsid`, `filter`
#'   - `af`, `ac`, `an`, `nhomalt` — global
#'   - `af_afr`, `af_amr`, `af_asj`, `af_eas`, `af_fin`, `af_mid`,
#'     `af_nfe`, `af_sas`, `af_remaining` — per-ancestry AFs
#'   - `grpmax`, `fafmax_faf95_max`, `fafmax_faf95_max_gen_anc`
#' @examples
#' \dontrun{
#' annotateGnomad(c(
#'   "chr1:45549863:A:ATATCTG",
#'   "chr1:45542866:TTTGAACTTCG:T"
#' ))
#' }
#' @export
annotateGnomad <- function(variant_ids,
                           dataset = c("genomes", "exomes"),
                           version = "4.1",
                           endpoint = "https://gnomad-public-us-east-1.s3.amazonaws.com",
                           cache_dir = NULL,
                           region_pad_bp = 0L) {
  dataset <- match.arg(dataset)
  stopifnot(
    is.character(variant_ids), length(variant_ids) >= 1L,
    is.character(endpoint), length(endpoint) == 1L,
    is.numeric(region_pad_bp), length(region_pad_bp) == 1L, region_pad_bp >= 0
  )
  parts <- .parse_variant_ids(variant_ids)

  # Group by chromosome and pull each range.
  groups <- split(seq_along(variant_ids), parts$chr)
  records <- vector("list", length(groups))
  names(records) <- names(groups)

  for (chr in names(groups)) {
    idx <- groups[[chr]]
    start <- max(1L, min(parts$pos[idx]) - as.integer(region_pad_bp))
    end   <- max(parts$pos[idx]) + as.integer(region_pad_bp)
    url <- .gnomad_vcf_url(endpoint = endpoint, dataset = dataset,
                           version = version, chr = chr)
    cli::cli_inform(c("i" = "Fetching {.val {chr}}:{start}-{end} from {.val {dataset}} v{version}."))
    rows <- .pull_region(url = url, chr = chr, start = start, end = end,
                         cache_dir = cache_dir, version = version,
                         dataset = dataset)
    records[[chr]] <- rows
  }

  all_rows <- do.call(rbind, records)
  if (is.null(all_rows) || nrow(all_rows) == 0L)
    all_rows <- .empty_gnomad_df()

  # Match each input variant to its gnomAD row (by full chr:pos:ref:alt).
  key_in  <- sprintf("%s:%d:%s:%s", parts$chr, parts$pos, parts$ref, parts$alt)
  key_out <- sprintf("%s:%d:%s:%s",
                     all_rows$chr, all_rows$pos, all_rows$ref, all_rows$alt)
  m <- match(key_in, key_out)

  tibble::tibble(
    variant_id                 = variant_ids,
    chr                        = parts$chr,
    pos                        = parts$pos,
    ref                        = parts$ref,
    alt                        = parts$alt,
    rsid                       = all_rows$rsid[m],
    filter                     = all_rows$filter[m],
    af                         = all_rows$af[m],
    ac                         = all_rows$ac[m],
    an                         = all_rows$an[m],
    nhomalt                    = all_rows$nhomalt[m],
    af_afr                     = all_rows$af_afr[m],
    af_amr                     = all_rows$af_amr[m],
    af_asj                     = all_rows$af_asj[m],
    af_eas                     = all_rows$af_eas[m],
    af_fin                     = all_rows$af_fin[m],
    af_mid                     = all_rows$af_mid[m],
    af_nfe                     = all_rows$af_nfe[m],
    af_sas                     = all_rows$af_sas[m],
    af_remaining               = all_rows$af_remaining[m],
    grpmax                     = all_rows$grpmax[m],
    fafmax_faf95_max           = all_rows$fafmax_faf95_max[m],
    fafmax_faf95_max_gen_anc   = all_rows$fafmax_faf95_max_gen_anc[m]
  )
}

# -- internals ---------------------------------------------------------

.gnomad_vcf_url <- function(endpoint, dataset, version, chr) {
  # e.g. .../release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz
  chr_str <- if (startsWith(chr, "chr")) chr else paste0("chr", chr)
  sprintf(
    "%s/release/%s/vcf/%s/gnomad.%s.v%s.sites.%s.vcf.bgz",
    sub("/$", "", endpoint), version, dataset, dataset, version, chr_str
  )
}

.parse_variant_ids <- function(ids) {
  parts <- strsplit(ids, ":", fixed = TRUE)
  if (any(lengths(parts) != 4L))
    cli::cli_abort(
      "All variant_ids must be {.val chr:pos:ref:alt}; bad: {.val {ids[lengths(parts) != 4L][1]}}"
    )
  m <- do.call(rbind, parts)
  chr <- ifelse(startsWith(m[, 1], "chr"), m[, 1], paste0("chr", m[, 1]))
  list(chr = chr, pos = as.integer(m[, 2]), ref = m[, 3], alt = m[, 4])
}

# Tab-separated VCF record parser: returns a data.frame with the
# columns we care about. Cheap, purpose-built — avoids pulling in
# VariantAnnotation just for AF extraction.
.parse_vcf_records <- function(lines) {
  if (length(lines) == 0L) return(.empty_gnomad_df())
  f <- do.call(rbind, strsplit(lines, "\t", fixed = TRUE))
  info <- f[, 8]

  get_tag <- function(tag) {
    pat <- sprintf("(^|;)%s=([^;]*)", tag)
    m <- regmatches(info, regexec(pat, info))
    vapply(m, function(v) if (length(v) >= 3) v[3] else NA_character_, character(1))
  }
  num <- function(tag) suppressWarnings(as.numeric(get_tag(tag)))
  int <- function(tag) suppressWarnings(as.integer(get_tag(tag)))

  data.frame(
    chr = f[, 1],
    pos = as.integer(f[, 2]),
    rsid = ifelse(f[, 3] == ".", NA_character_, f[, 3]),
    ref = f[, 4],
    alt = f[, 5],
    filter = f[, 7],
    af = num("AF"),
    ac = int("AC"),
    an = int("AN"),
    nhomalt = int("nhomalt"),
    af_afr = num("AF_afr"),
    af_amr = num("AF_amr"),
    af_asj = num("AF_asj"),
    af_eas = num("AF_eas"),
    af_fin = num("AF_fin"),
    af_mid = num("AF_mid"),
    af_nfe = num("AF_nfe"),
    af_sas = num("AF_sas"),
    af_remaining = num("AF_remaining"),
    grpmax = get_tag("grpmax"),
    fafmax_faf95_max = num("fafmax_faf95_max"),
    fafmax_faf95_max_gen_anc = get_tag("fafmax_faf95_max_gen_anc"),
    stringsAsFactors = FALSE
  )
}

.empty_gnomad_df <- function() {
  .parse_vcf_records(character())
  # ^ does not actually return an empty df because of do.call(rbind); build manually:
  data.frame(
    chr = character(), pos = integer(), rsid = character(), ref = character(),
    alt = character(), filter = character(), af = numeric(), ac = integer(),
    an = integer(), nhomalt = integer(),
    af_afr = numeric(), af_amr = numeric(), af_asj = numeric(), af_eas = numeric(),
    af_fin = numeric(), af_mid = numeric(), af_nfe = numeric(), af_sas = numeric(),
    af_remaining = numeric(), grpmax = character(),
    fafmax_faf95_max = numeric(), fafmax_faf95_max_gen_anc = character(),
    stringsAsFactors = FALSE
  )
}

.pull_region <- function(url, chr, start, end, cache_dir, version, dataset) {
  cache_path <- NULL
  if (!is.null(cache_dir)) {
    if (!dir.exists(cache_dir))
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    key <- sprintf("%s_%s_%s_%d_%d.rds", dataset, version, chr, start, end)
    cache_path <- file.path(cache_dir, key)
    if (file.exists(cache_path)) {
      cli::cli_inform(c("v" = "Cache hit: {.path {cache_path}}"))
      return(readRDS(cache_path))
    }
  }

  if (!requireNamespace("Rsamtools", quietly = TRUE))
    cli::cli_abort("Rsamtools is required for gnomAD annotation.")

  tf <- Rsamtools::TabixFile(url)
  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  raw <- Rsamtools::scanTabix(tf, param = gr)
  rows <- .parse_vcf_records(raw[[1]])

  if (!is.null(cache_path)) saveRDS(rows, cache_path)
  rows
}
