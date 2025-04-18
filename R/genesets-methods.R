genesets <- function(gs, ..., organism = NULL, strategy = NULL,
                     cache = NULL, verbose = TRUE) {
    UseMethod("genesets")
}

map_genesets <- function(gs, key_source, key_target,
                         annodb = NULL, organism = NULL) {
    assert_s3_class(gs, "enricher_genesets")
    assert_string(key_source, allow_empty = FALSE)
    assert_string(key_target, allow_empty = FALSE)
    if (vec_size(gs) == 0L) return(gs) # styler: off

    # Infer the annodb ----------------------------
    if (vec_unique_count(attr(gs, "organism")) > 1L) {
        cli::cli_abort(c(
            "Mapping genesets from multiple {.arg organism} values is not supported.",
            i = "Please map genesets for each organism separately before merging."
        ))
    }
    if (is.null(annodb)) {
        if (is.null(organism)) {
            organism <- vec_unique(attr(gs, "organism"))
            if (is.null(organism)) {
                cli::cli_abort(paste(
                    "Please specify {.arg annodb} to map {.arg gs} into",
                    "{.field {key_target}}, or specify {.arg organism}",
                    "so we can infer {.arg annodb} for you."
                ))
            }
        } else if (!.rlang_check_string(organism, allow_empty = FALSE)) {
            cli::cli_abort("{.arg organism} must be a single string")
        }

        annodb <- switch(organism,
            hsa = "org.Hs.eg.db",
            cli::cli_abort(c(
                "Cannot infer {.arg annodb} for {.field {organism}}",
                i = "Please provide {.arg annodb} manually"
            ))
        )
    }

    # Infer the annodb ----------------------------
    if (is.character(annodb)) {
        check_bioc_installed(annodb)
        annodb <- getExportedValue(annodb, annodb)
    }

    # mapping the the genes in genesets into keytype
    genesets_lapply(gs, function(geneset) {
        if (length(geneset) == 0L) return(geneset) # styler: off
        out <- AnnotationDbi::mapIds(
            annodb,
            keys = geneset,
            column = key_target,
            keytype = key_source
        )
        out[is.na(out) | out == ""] <- NA_character_
        out
    })
}

genesets_lapply <- function(gs, ...) {
    data <- vec_proxy(data)
    data$genesets <- lapply(data$genesets, ...)
    vec_restore(data, gs)
}

#' @export
genesets.character <- function(gs, ...) {
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
