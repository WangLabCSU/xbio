#' @importFrom rlang caller_call arg_match0
db <- function(download, ..., strategy, output, odir, ofile, description,
               verbose, arg_strategy = caller_arg(strategy),
               call = caller_call()) {
    if (!is.null(strategy)) {
        strategy <- arg_match0(
            strategy,
            c("read", "save", "download"),
            arg_nm = arg_strategy,
            error_call = call
        )
    }
    if (is.null(output)) {
        output <- file_path(odir, ofile, ext = "rds")
        if (identical(strategy, "save")) dir_create(odir, recursive = TRUE)
    } else {
        output <- paste0(output, ".rds")
    }
    if (is.null(strategy) || identical(strategy, "read")) {
        if (file.exists(output)) {
            out <- readRDS(output)
            if (isTRUE(verbose)) {
                cli::cli_inform(c(
                    "Using saved file: {.path {output}}",
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
                    "No file {.path {output}} found",
                    i = sprintf("Please save %s first", description)
                ),
                call = call
            )
        }
        strategy <- "save"
    }
    if (isTRUE(verbose)) cli::cli_inform(sprintf("Downloading %s", description))
    out <- download(...)
    if (identical(strategy, "save")) {
        if (isTRUE(verbose)) {
            cli::cli_inform(sprintf("Saving %s", description))
        }
        saveRDS(list(data = out, datetime = Sys.time()), file = output)
    }
    out
}

dbdir <- function(db) file_path(data_dir(), "db", db)
