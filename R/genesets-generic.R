#' Create Gene Sets
#'
#' @param gs Specify the gene sets to use. Accepted formats:
#' - `character`: A string indicating the source of gene sets. Acceptable values
#'   include `r oxford_or(c("go/GO", "kegg/KEGG"))`, which will download gene
#'   sets from the `GO` or `KEGG` databases. You may also specify the name of an
#'   [OrgDb](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb)
#'   package, which will call the `OrgDb` method. Alternatively, provide a file
#'   path to a `.gmt` or `.gmx` file to import custom gene sets.
#' - `data.frame`: A data frame with `2`, `3`, or `4` columns, representing gene
#'   set IDs, optional terms (if `4` columns), optional descriptions (if
#'   `3`/`4` columns), and associated genes or features.
#' - `list`: A list where each element is a character vector representing a gene
#'   set. A `description` attribute can be attached to each element, or a
#'   `descriptions` attribute to the entire list.
#' - `OrgDb`: An
#'   [OrgDb](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb)
#'   object used to retrieve GO gene sets.
#' @param ... Additional arguments passed to the method.
#' @inheritParams keggdb
#' @return A `enricher_genesets` object.
#' @export
genesets <- function(gs, ...) UseMethod("genesets")

#' @param select A character or integer vector select the columns to be used for
#' a data.frame input.
#' @export
#' @rdname genesets
genesets.data.frame <- function(gs, ..., select = NULL) {
    rlang::check_dots_empty()
    if (!is.null(select)) {
        select <- vec_as_location(
            select,
            n = ncol(gs),
            names = names(gs),
            missing = "error"
        )
        if (length(select) != 2L &&
            length(select) != 3L &&
            length(select) != 4L) {
            cli::cli_abort("{.arg select} must be of length 2, 3, or 4.")
        }
        gs <- gs[select]
    }
    vec_cast(gs, new_genesets())
}

#' @export
genesets.data.table <- function(gs, ...) genesets(as.data.frame(gs), ...)

#' @param ids A character vector of gene set IDs for list input.
#' @param terms A character vector of gene set terms for list input.
#' @param descriptions A character vector of gene set descriptions for a list
#' input.
#' @export
#' @rdname genesets
genesets.list <- function(gs, ..., ids = NULL, terms = NULL,
                          descriptions = NULL) {
    rlang::check_dots_empty()
    if (!is.null(ids)) {
        ids <- vec_cast(ids, character(), x_arg = "names")
        if (vec_any_missing(ids) || any(ids == "")) {
            cli::cli_abort("{.arg ids} cannot be missing or empty.")
        }
        names(gs) <- ids
    }
    if (!is.null(terms)) {
        terms <- vec_cast(terms, character())
        if (vec_size(gs) != vec_size(terms)) {
            cli::cli_abort(paste(
                "{.arg terms} must be",
                "the same length of the input list {.arg gs}."
            ))
        }
        attr(gs, "terms") <- terms
    }
    if (!is.null(descriptions)) {
        descriptions <- vec_cast(descriptions, character())
        if (vec_size(gs) != vec_size(descriptions)) {
            cli::cli_abort(paste(
                "{.arg descriptions} must be",
                "the same length of the input list {.arg gs}."
            ))
        }
        attr(gs, "descriptions") <- descriptions
    }
    vec_cast(gs, new_genesets())
}

#' @export
#' @rdname genesets
genesets.character <- function(gs, ..., verbose = TRUE) {
    if (!.rlang_check_string(gs, allow_empty = FALSE)) {
        cli::cli_abort("{.arg gs} must be a single non-empty string")
    }

    # File based method
    if (endsWith(tolower(gs), ".gmt")) {
        rlang::check_dots_empty()
        if (isTRUE(verbose)) {
            cli::cli_inform("Reading {.cls genesets} from {.path {gs}}")
        }
        return(read_gmt(path = gs, ...))
    } else if (endsWith(tolower(gs), ".gmx")) {
        rlang::check_dots_empty()
        if (isTRUE(verbose)) {
            cli::cli_inform("Reading {.cls genesets} from {.path {gs}}")
        }
        return(read_gmx(path = gs, ...))
    }

    # orgdb based method
    if (startsWith(gs, "org.") && endsWith(gs, ".db")) {
        if (isTRUE(verbose)) {
            cli::cli_inform(
                "Reading GO {.cls genesets} from {.pkg {gs}} package"
            )
        }
        check_bioc_installed(gs)
        return(genesets(gs = getExportedValue(gs, gs), ...))
    }

    # Network based method
    out <- switch(gs,
        GO = ,
        go = genesets_go(..., verbose = verbose),
        KEGG = ,
        kegg = genesets_kegg(..., verbose = verbose),
        # msigdbr = genesets_msigdbr(..., verbose = verbose),
        cli::cli_abort("No genesets for {gs}")
    )
    genesets(out)
}

#' @param keytype A single character string indicating which type of keys to
#' use.
#' @param ontology A character vector specifying the GO ontology to use. Must be
#' one or more of `r oxford_and(GO_ONTOLOGY)`.
#' @export
#' @rdname genesets
genesets.OrgDb <- function(gs, ..., keytype = "SYMBOL", ontology = NULL) {
    rlang::check_dots_empty()
    if (!is.null(ontology) && !all(ontology %in% GO_ONTOLOGY)) {
        cli::cli_abort(sprintf(
            "{.arg ontology} must be in {.field %s}",
            oxford_and(GO_ONTOLOGY, code = FALSE)
        ))
    }
    out <- AnnotationDbi::select(
        gs,
        keys = AnnotationDbi::keys(gs, keytype),
        columns = c("ONTOLOGYALL", "GOALL"),
        keytype = keytype
    )
    out <- vec_slice(out, !is.na(out$GO))
    if (!is.null(ontology)) {
        out <- vec_slice(out, out$ONTOLOGYALL %in% ontology)
    }
    genesets(out, select = c("GOALL", "ONTOLOGYALL", keytype))
}

#' @export
genesets.enricher_kegg_genesets <- function(gs, ...) {
    vec_cast(gs, new_genesets())
}

genesets_kegg <- function(organism = NULL, strategy = NULL, cache = NULL,
                          verbose = TRUE, ...) {
    rlang::check_dots_empty()
    keggdb(
        database = "genesets", organism = organism,
        strategy = strategy, cache = cache, verbose = verbose
    )
}

genesets_go <- function(organism = NULL, strategy = NULL, cache = NULL,
                        verbose = TRUE, ...) {
    rlang::check_dots_empty()
    godb(
        database = "genesets", organism = organism,
        strategy = strategy, cache = cache, verbose = verbose
    )
}
