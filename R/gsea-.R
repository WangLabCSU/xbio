gsea <- function(gs, object, method, ...,
                 min_size = NULL, max_size = NULL) {
    assert_s3_class(gs, "genesets")
    assert_number_whole(min_size, min = 1, allow_null = TRUE)
    assert_number_whole(min_size, min = 1, allow_null = TRUE)
    if (!all(keep <- list_sizes(gs) > 0L)) {
        cli::cli_warn("Removing {sum(!keep)} empty gene set{?s}")
        gs <- gs[keep]
    }
    gs <- genesets_lapply(gs, function(geneset) {
        geneset[!is.na(geneset) & geneset != ""]
    })
    if (!all(keep <- list_sizes(gs) > 0L)) {
        cli::cli_warn(paste(
            "Removing {sum(!keep)} invalid gene set{?s}",
            "(all are empty string or missing value)"
        ))
        gs <- gs[keep]
    }
    gs <- filter_genesets(
        gs,
        min_size = min_size,
        max_size = max_size
    )
    if (vec_size(gs) == 0L) cli::cli_abort("No gene sets {.arg gs} to use")
    gsea0(object, method, ..., gs = gs)
}

filter_genesets <- function(gs, min_size, max_size) {
    if (is.null(min_size) && is.null(max_size)) {
        return(gs)
    }
    sizes <- list_sizes(gs)
    if (!is.null(min_size) && !is.null(max_size)) {
        keep <- sizes >= min_size & sizes <= max_size
        out_pattern <- sprintf("[%d, %d]", min_size, max_size)
    } else if (!is.null(min_size)) {
        keep <- sizes >= min_size
        out_pattern <- sprintf("[%d, Inf)", min_size)
    } else {
        keep <- sizes <= max_size
        out_pattern <- sprintf("(0, %d]", max_size)
    }
    if (!all(keep)) {
        cli::cli_warn(sprintf(
            "Removing {sum(!keep)} gene set{?s} out of %s",
            out_pattern
        ))
        gs <- gs[keep]
    }
    gs
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
