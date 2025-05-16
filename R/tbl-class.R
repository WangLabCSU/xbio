new_tbl <- function(x, ..., class = NULL) {
    new_data_frame(x, ..., class = c(class, "tbl"))
}

new_gsea_result <- function(x, ...) new_tbl(x, ..., class = "xbio_gsea_result")

tbl_sum.xbio_gsea_result <- function(x, ...) {
    c("GSEA Result" = pillar::dim_desc(x))
}
