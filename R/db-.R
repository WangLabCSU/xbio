#' @importFrom rlang caller_call arg_match0
db <- function(download, ..., strategy, cache, cachedir, cachefile, description,
               verbose, arg_strategy = caller_arg(strategy),
               call = caller_call()) {
    if (!is.null(strategy)) {
        strategy <- arg_match0(
            strategy,
            c("read", "cache", "download"),
            arg_nm = arg_strategy,
            error_call = call
        )
    }
    if (!is.null(cache)) {
        cachedir <- getwd()
        cachefile <- cache
    }
    cached <- file_path(cachedir, cachefile, ext = "rds")
    if (is.null(strategy) || identical(strategy, "read")) {
        if (file.exists(cached)) {
            out <- readRDS(cached)
            if (isTRUE(verbose)) {
                cli::cli_inform(c(
                    "Using cached file: {.path {cached}}",
                    " " = sprintf(
                        "Snapshot date: %s",
                        as.character(.subset2(out, "datetime"), digits = 0L)
                    )
                ))
            }
            return(.subset2(out, "data"))
        } else if (identical(strategy, "read")) {
            cli::cli_abort(
                c(
                    "No cached file {.path {cached}} found",
                    i = sprintf("Please cache %s first", description)
                ),
                call = call
            )
        }
        strategy <- "cache"
    }
    if (isTRUE(verbose)) cli::cli_inform(sprintf("Downloading %s", description))
    out <- download(...)
    if (identical(strategy, "cache")) {
        dir_create(cachedir, recursive = TRUE)
        if (isTRUE(verbose)) {
            cli::cli_inform(sprintf("Caching %s", description))
        }
        saveRDS(list(data = out, datetime = Sys.time()), file = cached)
    }
    out
}

dbdir <- function(db) file_path(data_dir(), "db", db)
