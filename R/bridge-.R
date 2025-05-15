#' Bridge Functional Links Between Biological Processes
#'
#' `xbio()` / `bridge()` evaluate the functional association between source and
#' target biological process representations using various enrichment-based
#' methods. The two functions are equivalent in functionality but differ in the
#' ordering of their arguments.
#'
#' @param source The input source to be scored, depending on the selected
#' `method`.
#'  - **Over-Representation Analysis**: [`repr_threshold()`].
#'  - **GSEA gene-permutation algorithm**: [`repr_metrics()`].
#'  - **GSEA sample-permutation algorithm**: Accepts a matrix-like object.
#'    Although this method inherently relies on a ranked metric, it internally
#'    computes the ranking metric and performs permutations based on sample
#'    labels. Metric-related parameters should be specified within the `method`
#'    argument.
#'
#'    `Note`: While matrices are not conceptually equivalent to biological
#'    process representations, we provide [`repr_matrix()`] for stylistic
#'    consistency and user convenience.
#'  - **Rank-Rank Hypergeometric Overlap test**: [`repr_metrics()`].
#'
#' @param target The input target used for associating, depending on the
#' selected `method`.
#'  - **Over-Representation Analysis**: [`repr_genesets()`].
#'  - **GSEA gene-permutation algorithm**: [`repr_genesets()`].
#'  - **GSEA sample-permutation algorithm**: [`repr_genesets()`].
#'  - **Rank-Rank Hypergeometric Overlap test**: [`repr_metrics()`].
#'
#' @param method The scoring method to use:
#'  - **Over-Representation Analysis**: [ORA()]
#'  - **GSEA gene-permutation algorithm**: [`GSEAGene()`], [`GSEASimple()`],
#'    [`GSEAMultilevel()`], and [`GSEABroadGene()`].
#'  - **GSEA sample-permutation algorithm**: [`GSEABroad()`] and [GSEASample()].
#'  - **Rank-Rank Hypergeometric Overlap test**: [`RRHO()`].
#'
#' @details
#' In biological research, scoring methods are commonly used to evaluate the
#' association between two sets of biological processes. These may include
#' well-characterized processes (the *target*), poorly understood ones (the
#' *source*), or even two processes of unknown function. By quantifying their
#' association, we can infer potential biological relevance or shared
#' functionality. A strong association suggests functional similarity or
#' coordinated regulation.
#'
#' Biological processes are typically represented by one of the following:
#' - One or more gene sets ([`repr_genesets()`])
#' - A ranked metric ([`repr_metrics()`])
#' - A gene set based on a hard threshold ([`repr_threshold()`])
#'
#' @export
xbio <- function(source, target, method) bridge(method, source, target)

#' @export
#' @rdname xbio
bridge <- new_generic(
    "bridge", "method",
    function(method, source, target) {
        source <- repr_source(method, source)
        target <- repr_target(method, target)
        S7_dispatch()
    }
)

method(bridge, class_any) <- function(method, source, target) {
    cli::cli_abort(
        "No {.field bridge} method for {.obj_type_friendly {method}}"
    )
}

repr_source <- new_generic(
    "repr_source", c("method", "source"),
    function(method, source) S7_dispatch()
)

method(repr_source, list(class_missing, class_any)) <- function(method, source) {
    cli::cli_abort("{.arg method} must be provided")
}

method(repr_source, list(class_any, class_any)) <- function(method, source) {
    cli::cli_abort("Invalid {.arg method} provided")
}

repr_target <- new_generic(
    "repr_target", c("method", "target"),
    function(method, target) S7_dispatch()
)

method(repr_target, list(class_missing, class_any)) <- function(method, target) {
    cli::cli_abort("{.arg method} must be provided")
}

method(repr_target, list(class_any, class_any)) <- function(method, target) {
    cli::cli_abort("Invalid {.arg method} provided")
}
