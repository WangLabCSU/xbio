#' Get database from KEGG
#'
#' @param database A KEGG database (list available via
#' [`listDatabases()`][KEGGREST::listDatabases]). In addition, `"organism"` can
#' be used to retrieve all available organisms in the KEGG database, and
#' `"genesets"` can be used to obtain gene sets for functional enrichment
#' analysis.
#' @param organism A KEGG database (list available via
#' `keggdb("organism")$organism`.
#' @param strategy Character string specifying how to access the database and
#' whether to use caching:
#'  - `"read"`: Read from the saved file only. Error if not exist.
#'  - `"download"`: Download fresh data only. Do not save it.
#'  - `"save"`: Download data and save it locally.
#'
#' By default, `NULL` is used, which means read from saved file if available;
#' otherwise, download and save the data.
#'
#' @param save A string specifying the file path to read from or save to.
#' @param verbose A single logical value indicates whether the process should be
#' verbose.
#' @return A list of contents for each entries.
#' @importFrom rlang arg_match0
#' @export
keggdb <- function(database, organism = NULL, strategy = NULL, save = NULL,
                   verbose = TRUE) {
    check_bioc_installed("KEGGREST", "to download from KEGG")
    database <- arg_match0(database, c(
        "organism", "genesets", KEGGREST::listDatabases()
    ))
    if (identical(database, "organism")) {
        if (!is.null(organism)) {
            cli::cli_abort(paste(
                "{.arg organism} cannot be used for",
                "{.field organism} database"
            ))
        }
        description <- "{.field organism} database in KEGG"
        odir <- dbdir("KEGG")
        ofile <- "organism"
    } else {
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
        odir <- file_path(dbdir("KEGG"), database)
        ofile <- organism
    }
    db(
        keggdb_download, database, organism,
        strategy = strategy, output = save,
        odir = odir, ofile = ofile,
        description = description, verbose = verbose
    )
}

keggdb_download <- function(database, organism) {
    if (identical(database, "organism")) {
        out <- KEGGREST::keggList("organism")
        # https://stackoverflow.com/questions/6819804/convert-a-matrix-to-a-list-of-column-vectors
        # For performances
        out <- lapply(seq_len(ncol(out)), function(i) out[, i, drop = TRUE])
        structure(out, class = sprintf("%s_kegg_%s", pkg_nm(), database))
    } else if (identical(database, "genesets")) {
        pathway2genes <- KEGGREST::keggLink(organism, "pathway")
        pathways <- sub("^[^:]+:", "", names(pathway2genes))
        terms <- KEGGREST::keggList("pathway", organism)
        genes <- sub("^[^:]+:", "", pathway2genes)
        new_data_frame(list(
            ids = pathways,
            terms = terms[pathways],
            genes = genes
        ), class = sprintf("%s_kegg_genesets", pkg_nm()))
    } else {
        items <- KEGGREST::keggList(database, organism)
        out <- keggdb_get0(names(items))
        structure(out, class = sprintf("%s_kegg_%s", pkg_nm(), database))
    }
}

methods::setOldClass("xbio_kegg_genesets")

S3_kegg_genesets <- new_S3_class("xbio_kegg_genesets")

#' Finds entries in a given database
#'
#' @param ... One or more keywords, or a range of integers representing
#' molecular weights. If identifiers included are not known to KEGG, the
#' results will not contain any information about those identifiers.
#' @inheritParams KEGGREST::keggFind
#' @importFrom rlang arg_match0
keggdb_find <- function(database, ..., option = NULL) {
    check_bioc_installed("KEGGREST", "to query from KEGG")
    database <- arg_match0(database, KEGGREST::listDatabases())
    query <- as.character(unlist(rlang::list2(...), use.names = FALSE))
    if (is.null(option)) {
        KEGGREST::keggFind(database = database, query = query)
    } else {
        KEGGREST::keggFind(database = database, query = query, option = option)
    }
}

#' Retrieves given database entries
#'
#' @param ... One or more KEGG identifiers.
#' @inheritParams KEGGREST::keggGet
keggdb_get <- function(..., option = NULL) {
    check_bioc_installed("KEGGREST", "to query from KEGG")
    query <- as.character(unlist(rlang::list2(...), use.names = FALSE))
    keggdb_get0(query, option = option)
}

keggdb_get0 <- function(query, option = NULL) {
    out <- vector("list", length(query))
    names(out) <- query
    index <- seq_along(out)
    groups <- split(index, (index %/% 10) + 1L)
    progress_id <- cli::cli_progress_bar(
        name = "keggGet",
        format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} [{cli::pb_rate}] | {cli::pb_eta_str}",
        format_done = "Get from KEGG for {.val {cli::pb_total}} quer{?y/ies} in {cli::pb_elapsed}",
        total = length(out),
        clear = FALSE
    )

    for (idx in groups) {
        if (is.null(option)) {
            out[idx] <- KEGGREST::keggGet(query[idx])
        } else {
            out[idx] <- KEGGREST::keggGet(query[idx], option = option)
        }
        cli::cli_progress_update(inc = length(idx), id = progress_id)
    }
    out
}
