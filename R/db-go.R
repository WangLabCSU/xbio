godb <- function(database, organism = NULL, strategy = NULL, cache = NULL,
                 verbose = TRUE) {
    organism <- check_organism(organism)
    if (identical(database, "genesets")) {
        description <- sprintf(
            "{.field %s} for {.field %s} organism in KEGG",
            database, organism
        )
    } else {
        description <- sprintf(
            "{.field %s} database for {.field %s} organism in KEGG",
            database, organism
        )
    }
    cachedir <- file_path(dbdir("KEGG"), database)
    cachefile <- organism
    db(
        godb_download, database, organism,
        strategy = strategy, cache = cache,
        cachedir = cachedir, cachefile = cachefile,
        description = description, verbose = verbose
    )
}

godb_download <- function(database, organism) {
    if (identical(database, "organism")) {
        out <- KEGGREST::keggList("organism")
        # https://stackoverflow.com/questions/6819804/convert-a-matrix-to-a-list-of-column-vectors
        # For performances
        out <- lapply(seq_len(ncol(out)), function(i) out[, i, drop = TRUE])
        structure(out, class = sprintf("%s_kegg_%s", pkg_nm(), database))
    } else if (identical(database, "genesets")) {
        pathway2genes <- KEGGREST::keggLink(organism, "pathway")
        pathways <- sub("^[^:]+:", "", names(pathway2genes))
        descriptions <- KEGGREST::keggList("pathway", organism)
        genes <- sub("^[^:]+:", "", pathway2genes)
        new_data_frame(list(
            pathways = pathways,
            descriptions = descriptions[pathways],
            genes = genes
        ), class = sprintf("%s_kegg_genesets", pkg_nm()))
    } else {
        items <- KEGGREST::keggList(database, organism)
        out <- keggdb_get0(names(items))
        structure(out, class = sprintf("%s_kegg_%s", pkg_nm(), database))
    }
}
