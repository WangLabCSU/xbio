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

check_remote_installed <- function(pkg, remote, reason = NULL, ...) {
    rlang::check_installed(
        pkg,
        reason = reason,
        ...,
        action = function(pkgs, ...) {
            if (is_installed("pak")) {
                getExportedValue("pak", "pkg_install")(remote, ask = FALSE, ...)
            } else if (is_installed("remotes")) {
                getExportedValue("remotes", "install_github")(remote, ...)
            } else {
                choosed <- utils::menu(
                    c("pak", "remotes"),
                    title = paste(
                        "Would you like to install `pak`/`remotes`",
                        "in order to install", oxford_and(remote)
                    )
                )
                if (choosed == 1L) {
                    utils::install.packages("pak")
                    getExportedValue("pak", "pkg_install")(
                        remote, ask = FALSE, ...
                    )
                } else if (choosed == 2L) {
                    utils::install.packages("remotes")
                    getExportedValue("remotes", "install_github")(remote, ...)
                } else {
                    invokeRestart("abort")
                }
            }
        }
    )
}

RUST_CALL <- .Call

#' @keywords internal
rust_method <- function(class, method, ...) {
    rust_call(sprintf("%s__%s", class, method), ...)
}

#' @keywords internal
rust_call <- function(.NAME, ...) {
    # call the function
    out <- RUST_CALL(sprintf("wrap__%s", .NAME), ...)

    # propagate error from rust --------------------
    if (!inherits(out, "extendr_result")) return(out) # styler: off
    if (is.null(.subset2(out, "ok"))) cli::cli_abort(.subset2(out, "err"))
    .subset2(out, "ok")
}
