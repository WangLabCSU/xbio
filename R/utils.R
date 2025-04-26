as_data_frame <- function(x) {
    if (is_installed("tibble")) {
        getExportedValue("tibble", "as_tibble")(x, .name_repair = "minimal")
    } else {
        as.data.frame(
            x,
            row.names = FALSE,
            make.names = FALSE,
            stringsAsFactors = FALSE
        )
    }
}

#' @importFrom rlang caller_call
check_organism <- function(organism, call = caller_call()) {
    if (is.null(organism)) {
        organism <- "hsa"
    } else {
        assert_string(organism, allow_empty = FALSE, call = call)
        available_organisms <- keggdb("organism", verbose = FALSE)
        if (any(tcodes <- organism == .subset2(available_organisms, 1L))) {
            organism <- .subset2(available_organisms, 2L)[which(tcodes)]
        } else if (!any(organism == .subset2(available_organisms, 2L))) {
            cli::cli_abort(
                "Cannot found {.field {organism}} organism",
                call = call
            )
        }
    }
    organism
}

check_bioc_installed <- function(pkg, reason = NULL, ...) {
    rlang::check_installed(
        pkg,
        reason = reason,
        ...,
        action = function(pkgs, ...) {
            if (is_installed("pak")) {
                getExportedValue("pak", "pkg_install")(pkgs, ask = FALSE, ...)
            } else if (is_installed("BiocManager")) {
                getExportedValue("BiocManager", "install")(pkgs, ...)
            } else {
                choosed <- utils::menu(
                    c("pak", "BiocManager"),
                    title = paste(
                        "Would you like to install `pak`/`BiocManager`",
                        "in order to install", oxford_and(pkgs)
                    )
                )
                if (choosed == 1L) {
                    utils::install.packages("pak")
                    getExportedValue("pak", "pkg_install")(
                        pkgs, ask = FALSE, ...
                    )
                } else if (choosed == 2L) {
                    utils::install.packages("BiocManager")
                    getExportedValue("BiocManager", "install")(pkgs, ...)
                } else {
                    invokeRestart("abort")
                }
            }
        }
    )
}

#' @keywords internal
call_rust_method <- function(class, method, ...) {
    call_rust_fn(sprintf("wrap__%s__%s", class, method), ...)
}

#' @keywords internal
call_rust_fn <- function(.NAME, ...) {
    # call the function
    out <- .standalone_types_check_assert_call(.NAME, ...)

    # propagate error from rust --------------------
    rust_unwrap(out)
}

rust_unwrap <- function(out) {
    if (!inherits(out, "extendr_result")) return(out) # styler: off
    if (is.null(.subset2(out, "ok"))) stop(.subset2(out, "err"), call. = FALSE)
    .subset2(out, "ok")
}
