#' Get database from KEGG
#'
#' @param database A KEGG database (list available via
#' [`listDatabases()`][KEGGREST::listDatabases]). In addition, `"organism"` can
#' be used to get all available organisms in KEGG database.
#' @param organism A KEGG database (list available via
#' `keggdb("organism")$organism`.
#' @param strategy Character string specifying how to access the database and
#' whether to use caching:
#'  - `"read"`: Read from cache only. Error if not cached.
#'  - `"download"`: Download fresh data. Do not cache it.
#'  - `"cache"`: Download data and cache it locally.
#'
#' By default, `NULL` is used, which means read from cache if available;
#' otherwise, download and cache the data.
#'
#' @param verbose A single logical value indicates whether the process should be
#' verbose.
#' @return A list of contents for each entries.
#' @importFrom rlang arg_match0
#' @export
keggdb <- function(database, organism = NULL, strategy = NULL, verbose = TRUE) {
    rlang::check_installed("KEGGREST", "to download from KEGG")
    database <- arg_match0(
        database, c("organism", "gsea", KEGGREST::listDatabases())
    )
    if (identical(database, "organism")) {
        description <- "{.field organism} database in KEGG"
        if (!is.null(organism)) {
            cli::cli_abort(paste(
                "{.arg organism} cannot be used for", description
            ))
        }
        cachedir <- dbdir("KEGG")
        cachefile <- "organism"
    } else {
        organism <- kegg_check_organism(organism)
        if (identical(database, "gsea")) {
            description <- sprintf(
                "{.field %s} pathways for {.field %s} organism in KEGG",
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
    }
    db(
        keggdb_download, database, organism,
        cachedir = cachedir, cachefile = cachefile, strategy = strategy,
        description = description, verbose = verbose
    )
}

#' @importFrom rlang caller_call
kegg_check_organism <- function(organism, call = caller_call()) {
    if (is.null(organism)) {
        organism <- "hsa"
    } else {
        available_organisms <- keggdb("organism", verbose = FALSE)
        if (any(tcodes <- organism == .subset2(available_organisms, 1L))) {
            organism <- .subset2(available_organisms, 2L)[which(tcodes)]
        } else if (!any(organism == .subset2(available_organisms, 2L))) {
            cli::cli_abort(
                "Cannot found {.field {organism}} organism",
                call = call
            )
        }
    }
    organism
}

keggdb_download <- function(database, organism) {
    if (identical(database, "organism")) {
        out <- apply(
            KEGGREST::keggList("organism"), 2L, identity,
            simplify = FALSE
        )
        structure(out, class = sprintf("%s_kegg_%s", pkg_nm(), database))
    } else if (identical(database, "gsea")) {
        pathway2genes <- KEGGREST::keggLink(organism, "pathway")
        pathways <- sub("^[^:]+:", "", names(pathway2genes), fixed = TRUE)
        descriptions <- KEGGREST::keggList("pathway", organism)
        genes <- sub("^[^:]+:", "", pathway2genes, fixed = TRUE)
        new_data_frame(list(
            pathways = pathways,
            descriptions = descriptions[pathways],
            genes = genes
        ), class = sprintf("%s_kegg_gsea", pkg_nm()))
    } else {
        items <- KEGGREST::keggList(database, organism)
        out <- keggdb_get0(names(items))
        structure(out, class = sprintf("%s_kegg_%s", pkg_nm(), database))
    }
}

#' Finds entries in a given database
#'
#' @param ... One or more keywords, or a range of integers representing
#' molecular weights. If identifiers included are not known to KEGG, the
#' results will not contain any information about those identifiers.
#' @inheritParams KEGGREST::keggFind
keggdb_find <- function(database, ..., option = NULL) {
    rlang::check_installed("KEGGREST", "to query from KEGG")
    database <- rlang::arg_match0(database, KEGGREST::listDatabases())
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
    rlang::check_installed("KEGGREST", "to query from KEGG")
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
    if (is.null(option)) {
        for (idx in groups) {
            out[idx] <- KEGGREST::keggGet(query[idx])
            cli::cli_progress_update(inc = length(idx), id = progress_id)
        }
    } else {
        for (idx in groups) {
            out[idx] <- KEGGREST::keggGet(query[idx], option = option)
            cli::cli_progress_update(inc = length(idx), id = progress_id)
        }
    }
    out
}
