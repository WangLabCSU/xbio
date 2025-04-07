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
