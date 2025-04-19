gsea <- function(gs, object, method, ...,
                 min_size = NULL, max_size = NULL) {
    assert_s3_class(gs, "enricher_genesets")
    if (!all(keep <- list_sizes(gs) > 0L)) {
        cli::cli_warn("Removing {sum(!keep)} empty gene set{?s}")
        gs <- gs[keep]
    }
    gs <- gs_trim(gs)
    if (!all(keep <- list_sizes(gs) > 0L)) {
        cli::cli_warn(paste(
            "Removing {sum(!keep)} invalid gene set{?s}",
            "(all are empty string or missing value)"
        ))
        gs <- gs[keep]
    }
    gs <- gs_filter(
        gs,
        min_size = min_size,
        max_size = max_size
    )
    if (vec_size(gs) == 0L) cli::cli_abort("No gene sets to use")
    gsea0(object, method, ..., gs = gs)
}

gsea0 <- new_generic(
    "gsea0", c("object", "method"),
    function(object, method, ..., gs) S7_dispatch()
)

method(gsea0, list(class_any, class_any)) <- function(object, method, ..., gs) {
    cli::cli_abort(paste(
        "No gsea method for the combined signature:",
        "{.obj_type_friendly {object}} and {.obj_type_friendly {method}}"
    ))
}
