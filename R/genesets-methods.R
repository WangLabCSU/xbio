genesets <- function(gs, ...) UseMethod("genesets")

#' @export
genesets.character <- function(gs, ..., organism = NULL, strategy = NULL,
                               cache = NULL, verbose = TRUE) {
    genesets_db <- vec_unique(gs)
    out <- vector("list", length(genesets_db))
    names(out) <- genesets_db
    for (dbname in genesets_db) {
        out[[dbname]] <- switch(dbname,
            GO = genesets_go(...),
            KEGG = genesets_kegg(...),
            cli::cli_abort("NO genesets for {x}")
        )
    }
    unlist(out, recursive = FALSE, use.names = TRUE)
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
