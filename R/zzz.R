.onLoad <- function(...) {
    s3_register("pillar::tbl_sum", "xbio_gsea_result")
    methods_register()
}
