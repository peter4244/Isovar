#' Load user-level secrets from a KEY=VALUE file
#'
#' Reads a simple shell-style `KEY=VALUE` file (one per line; blank lines
#' and `#` comments ignored) and sets each pair as an environment
#' variable for the current R session via [Sys.setenv()]. Existing
#' environment variables are **not** overwritten by default, so an
#' explicit `Sys.setenv()` or shell-exported value takes precedence
#' over the file. Pass `overwrite = TRUE` to force reload.
#'
#' The canonical location is `~/.config/isovar/secrets.env` — see the
#' README there for what belongs in it and how to manage rotation.
#'
#' @param path Path to the secrets file. Default `~/.config/isovar/secrets.env`.
#' @param overwrite If `TRUE`, overwrite existing env vars. Default `FALSE`.
#' @param quiet If `TRUE`, suppress the "loaded N keys" message.
#' @return Invisibly, a character vector of keys that were set.
#' @export
loadIsovarSecrets <- function(path = "~/.config/isovar/secrets.env",
                              overwrite = FALSE,
                              quiet = FALSE) {
  path <- path.expand(path)
  if (!file.exists(path)) {
    if (!quiet)
      cli::cli_inform(c("i" = "No secrets file at {.path {path}}; nothing loaded."))
    return(invisible(character()))
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!grepl("^\\s*#", lines)]

  set_keys <- character()
  for (ln in lines) {
    m <- regmatches(ln, regexec("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.*?)\\s*$", ln))[[1]]
    if (length(m) < 3L) next
    key <- m[2]; val <- m[3]
    # Strip optional surrounding quotes.
    val <- sub('^"(.*)"$',  "\\1", val)
    val <- sub("^'(.*)'$",  "\\1", val)
    if (!overwrite && nzchar(Sys.getenv(key, unset = ""))) next
    args <- list(val); names(args) <- key
    do.call(Sys.setenv, args)
    set_keys <- c(set_keys, key)
  }

  if (!quiet) {
    if (length(set_keys))
      cli::cli_inform(c("v" = "Loaded {length(set_keys)} secret{?s} from {.path {path}}: {.envvar {set_keys}}"))
    else
      cli::cli_inform(c("i" = "No new secrets loaded from {.path {path}} (already set or file empty)."))
  }
  invisible(set_keys)
}
