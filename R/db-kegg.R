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
#' @param link A KEGG organism code (e.g., `"hsa"` for human). A full list is
#' available via `keggdb("organism")$organism`. This will link database using
#' gene symbols for the specified organism. Alternatively, a KEGG database name
#' can be provided to link pathway/disease associations based on database entry
#' IDs.
#'
#' the following database name can be linked to **pathway** database:
#' `r oxford_and(setdiff(KEGGREST::listDatabases(), KEGG_NO_PATHWAY))`.
#'
#' the following database name can be linked to **disease** database:
#' `r oxford_and(setdiff(KEGGREST::listDatabases(), KEGG_NO_DISEASE))`.
#'
#' @param database genesets database, a string of `"pathway"` or `"disease"`.
#' @inheritParams keggdb
#' @return A data frame
#' @export
kegg_genesets <- function(link = NULL, database = NULL, strategy = NULL,
                          save = NULL, verbose = TRUE) {
    assert_string(link, allow_empty = FALSE, allow_null = TRUE)
    link <- link %||% "hsa"
    if (is.null(database)) {
        database <- "pathway"
    } else {
        database <- rlang::arg_match0(database, c("disease", "pathway"))
    }
    if ((identical(database, "pathway") && any(link == KEGG_NO_PATHWAY)) ||
        (identical(database, "disease") && any(link == KEGG_NO_DISEASE))) {
        cli::cli_abort(sprintf(
            "Cannot link {.field %s} with KEGG {.field %s} database",
            link, database
        ))
    } else if (any(link == KEGGREST::listDatabases())) {
        organism <- NULL
        description <- sprintf("{.field %s}", link)
    } else {
        available_organisms <- keggdb("organism", verbose = FALSE)
        if (any(tcodes <- link == .subset2(available_organisms, 1L))) {
            link <- .subset2(available_organisms, 2L)[which(tcodes)]
        } else if (!any(link == .subset2(available_organisms, 2L))) {
            cli::cli_abort(sprintf(
                "Cannot find {.field %s} %s database",
                link, database
            ))
        }
        organism <- link
        if (identical(database, "disease") && !identical(organism, "hsa")) {
            cli::cli_abort(paste(
                "Only {.val hsa} {.arg organism} can be used for",
                "{.field disease} database"
            ))
        }
        description <- sprintf("{.field %s} organism", organism)
    }
    ofile <- sprintf("%s_%s", link, database)
    odir <- file_path(dbdir("KEGG"), "genesets")
    db(
        kegg_genesets_download,
        database = database, organism = organism, link = link,
        strategy = strategy, output = save,
        odir = odir, ofile = ofile,
        description = sprintf("KEGG {.field %s} for %s", database, description),
        verbose = verbose
    )
}

kegg_genesets_download <- function(database, organism, link) {
    ids2features <- KEGGREST::keggLink(link, database)
    if (identical(database, "disease")) {
        terms <- KEGGREST::keggList(database)
    } else if (is.null(organism)) { # for pathway database
        terms <- KEGGREST::keggList(database)
        if (any(link == c("ko", "reaction", "enzyme"))) {
            ids2features <- vec_slice(
                ids2features, grepl("^path:map", names(ids2features))
            )
        }
    } else {
        terms <- KEGGREST::keggList(database, organism)
    }
    ids <- sub("^[^:]+:", "", names(ids2features))
    features <- sub("^[^:]+:", "", ids2features)
    new_data_frame(
        list(
            ids = ids,
            terms = terms[ids],
            features = features
        ),
        class = sprintf("%s_kegg_genesets", pkg_nm())
    )
}

methods::setOldClass("xbio_kegg_genesets")

S3_kegg_genesets <- new_S3_class("xbio_kegg_genesets")

KEGG_NO_PATHWAY <- c(
    "pathway", "disease", "brite", "genome", "vg", "ag",
    "dgroup", "environ", "genes", "ligand", "kegg"
)

KEGG_NO_DISEASE <- c(
    "pathway", "disease", "genome", "module", "vg", "ag",
    "compound", "glycan", "reaction", "rclass", "enzyme",
    "dgroup", "environ", "genes", "ligand", "kegg"
)
