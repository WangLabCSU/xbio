#' Bridge Functional Links Between Biological Processes
#'
#' @description
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
#' @param source The input source to be scored, depending on the selected
#' `method`.
#'  - Over-Representation Analysis: [`repr_threshold()`].
#'  - GSEA gene-permutation algorithm: [`repr_metrics()`].
#'  - GSEA sample-permutation algorithm: A matrix-like object. Although this
#'    method also utilizes [`repr_metrics()`], it internally computes ranking
#'    metrics and performs permutations based on sample labels. Metric-related
#'    arguments should be specified within the `method` argument. 
#' 
#'    `Note`: While matrices are not conceptually equivalent to biological
#'    process representations, we provide [`repr_matrix()`] for stylistic
#'    consistency and user convenience.
#'  - Rank-Rank Hypergeometric Overlap test: [`repr_metrics()`].
#'
#' @param target The input target used for associating, depending on the
#' selected `method`.
#'  - Over-Representation Analysis: [`repr_genesets()`].
#'  - GSEA gene-permutation algorithm: [`repr_genesets()`].
#'  - GSEA sample-permutation algorithm: [`repr_genesets()`].
#'  - Rank-Rank Hypergeometric Overlap test: [`repr_metrics()`].
#'
#' @param method The scoring method to use:
#'  - Over-Representation Analysis: [ORA()]
#'  - GSEA gene-permutation algorithm: [`GSEAGene()`], [`GSEASimple()`],
#'    [`GSEAMultilevel()`], and [`GSEABroadGene()`].
#'  - GSEA sample-permutation algorithm: [`GSEABroad()`] and [GSEASample()].
#'  - Rank-Rank Hypergeometric Overlap test: [`RRHO()`].
#' @export
methods::setGeneric(
    "bridge",
    signature = "method",
    function(source, target, method) {
        source <- repr_source(method, source)
        target <- repr_target(method, target)
        standardGeneric("bridge")
    }
)

methods::setMethod("bridge", "ANY", function(source, target, method) {
    cli::cli_abort(
        "No {.field bridge} method for {.obj_type_friendly {method}}"
    )
})

repr_source <- new_generic(
    "repr_source", "method",
    function(method, source) S7_dispatch()
)

method(repr_source, class_missing) <- function(method, source) {
    cli::cli_abort("{.arg method} must be provided")
}


method(repr_source, class_any) <- function(method, source) {
    cli::cli_abort("Invalid {.arg method} provided")
}

repr_target <- new_generic(
    "repr_target", "method",
    function(method, target) S7_dispatch()
)

method(repr_target, class_missing) <- function(method, target) {
    cli::cli_abort("{.arg method} must be provided")
}

method(repr_target, class_any) <- function(method, target) {
    cli::cli_abort("Invalid {.arg method} provided")
}
