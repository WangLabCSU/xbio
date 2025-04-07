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
