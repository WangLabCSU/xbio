.onLoad <- function(...) {
    s3_register("pillar::tbl_sum", "xbio_gsea_result")
    s3_register("pillar::tbl_sum", "xbio_kegg_genesets")
    methods_register()
}
