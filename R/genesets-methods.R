genesets <- function(gs, ...) UseMethod("genesets")

#' @export
genesets.character <- function(gs, ...) {
    gs <- arg_match0(gs, c("go", "kegg"))
    out <- switch(gs,
        go = genesets_go(...),
        kegg = genesets_kegg(...),
        cli::cli_abort("NO genesets for {gs}")
    )
    genesets(out)
}

genesets_kegg <- function(organism = NULL, strategy = NULL, cache = NULL,
                          verbose = TRUE, ...) {
    keggdb(
        database = "genesets", organism = organism,
        strategy = strategy, cache = cache, verbose = verbose
    )
}

genesets_go <- function(organism = NULL, strategy = NULL, cache = NULL,
                        verbose = TRUE, ...) {
    godb(
        database = "genesets", organism = organism,
        strategy = strategy, cache = cache, verbose = verbose
    )
}

#' @export
genesets.data.frame <- function(gs, ...) vec_cast(gs, new_genesets())

#' @export
genesets.list <- function(gs, ...) vec_cast(gs, new_genesets())

#' @export
genesets.enricher_kegg_genesets <- function(gs, ...) {
    vec_cast(gs, new_genesets())
}
