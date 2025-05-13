#' Get database from KEGG
#'
#' @param database A KEGG database (list available via
#' [`listDatabases()`][KEGGREST::listDatabases]). In addition, `"organism"` can
#' be used to retrieve all available organisms in the KEGG database.
#' @param organism A KEGG organism code (list available via
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
#' @param save A string specifying the file path to read from or save to. File
#' extension `".rds"` will be automatically added.
#' @param verbose A single logical value indicates whether the process should be
#' verbose.
#' @return A list of contents for each entries.
#' @importFrom rlang arg_match0
#' @export
keggdb <- function(database, organism = NULL, strategy = NULL, save = NULL,
                   verbose = TRUE) {
    check_bioc_installed("KEGGREST", "to download from KEGG")
    database <- arg_match0(database, c("organism", KEGGREST::listDatabases()))
    if (any(database == c("organism", "ko"))) {
        if (!is.null(organism)) {
            cli::cli_abort(sprintf(
                "{.arg organism} cannot be used for {.field %s} database",
                database
            ))
        }
        description <- sprintf("{.field %s} database of KEGG", database)
        odir <- dbdir("KEGG")
        ofile <- database
    } else {
        organism <- check_organism(organism)
        description <- sprintf(
            "{.field %s} database for {.field %s} organism of KEGG",
            database, organism
        )
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
    } else {
        if (is.null(organism)) {
            out <- KEGGREST::keggList(database)
        } else {
            out <- KEGGREST::keggList(database, organism)
        }
        out <- keggdb_get0(names(out))
        structure(out, class = sprintf("%s_kegg_%s", pkg_nm(), database))
    }
}

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

#' Retrieve KEGG Gene Sets or Pathways
#'
#' Retrieves gene sets or pathways from the KEGG database for a given organism.
#'
#' @param database A KEGG organism code (e.g., `"hsa"` for human). A full list
#' is available via `keggdb("organism")$organism`. This will retrieve pathways
#' using gene symbols for the specified organism. Alternatively, a KEGG database
#' name can be provided to retrieve pathway associations based on database entry
#' IDs.
#'
#' **Note**: the following databases can be linked to pathways:
#' `r oxford_and(setdiff(KEGGREST::listDatabases(), KEGG_NO_PATHWAY))`.
#' @inheritParams keggdb
#' @return A data frame
#' @export
kegg_pathway <- function(database = "hsa", strategy = NULL,
                         save = NULL, verbose = TRUE) {
    assert_string(database, allow_empty = FALSE, allow_null = TRUE)
    database <- database %||% "hsa"
    if (any(database == KEGG_NO_PATHWAY)) {
        cli::cli_abort(sprintf(
            "Cannot link {.field %s} database with KEGG {.field pathway}",
            database
        ))
    } else if (any(database == KEGGREST::listDatabases())) {
        description <- sprintf("{.field %s}", database)
        ofile <- sprintf("%s_pathway", database)
        organism <- NULL
    } else {
        available_organisms <- keggdb("organism", verbose = FALSE)
        if (any(tcodes <- database == .subset2(available_organisms, 1L))) {
            database <- .subset2(available_organisms, 2L)[which(tcodes)]
        } else if (!any(database == .subset2(available_organisms, 2L))) {
            cli::cli_abort("Cannot found {.field {database}} pathway database")
        }
        organism <- database
        description <- sprintf("{.field %s} organism", organism)
        database <- NULL
        ofile <- sprintf("%s_pathway", organism)
    }
    odir <- file_path(dbdir("KEGG"), "genesets")
    db(
        kegg_pathway_download,
        database = database, organism = organism,
        strategy = strategy, output = save,
        odir = odir, ofile = ofile,
        description = sprintf("KEGG {.field pathway} for %s", description),
        verbose = verbose
    )
}

kegg_pathway_download <- function(database, organism) {
    if (is.null(database)) {
        pathway2features <- KEGGREST::keggLink(organism, "pathway")
        terms <- KEGGREST::keggList("pathway", organism)
    } else {
        pathway2features <- KEGGREST::keggLink(database, "pathway")
        terms <- KEGGREST::keggList("pathway")
        if (!any(database == c("rclass", "disease"))) {
            pathway2features <- vec_slice(
                pathway2features,
                grepl("^path:map", names(pathway2features))
            )
        }
    }
    pathways <- sub("^[^:]+:", "", names(pathway2features))
    features <- sub("^[^:]+:", "", pathway2features)
    new_data_frame(
        list(
            ids = pathways,
            terms = terms[pathways],
            features = features
        ),
        class = sprintf("%s_kegg_genesets", pkg_nm())
    )
}

methods::setOldClass("xbio_kegg_genesets")

S3_kegg_genesets <- new_S3_class("xbio_kegg_genesets")

KEGG_NO_PATHWAY <- c(
    "pathway", "brite", "genome", "vg", "ag",
    "dgroup", "environ", "genes", "ligand", "kegg"
)
