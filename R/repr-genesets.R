#' Create genesets representation
#'
#' @description
#' `repr_genesets()` is a wrapper around `genesets()` that additionally cleans
#' and filters the gene sets. It removes invalid entries, such as gene sets with
#' no genes or those containing only missing values. This function ensures the
#' returned object is clean and ready for downstream analysis.
#'
#' All analyses that depend on the gene set representation will use
#' `repr_genesets()` with its default arguments. If you want to customize
#' filtering, such as setting thresholds on gene set size (`min_size` /
#' `max_size`), you must call `repr_genesets()`/`gs_filter()` manually.
#'
#' @param gs Specify the gene sets to use. Accepted formats:
#' - `character`: A string indicating the source of gene sets. Valid values
#'   include:
#'   * `"go"` / `"GO"` or `"kegg"` / `"KEGG"`: Download gene sets from GO or
#'     KEGG using [`godb()`] or [`kegg_genesets()`]. Optional parameters include
#'     `link`, `database`, `strategy`, and `save`. See [`kegg_genesets()`] and
#'     [`godb()`] for details.
#'   * The name of an [`OrgDb`](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb) package, used to extract GO gene sets.
#'   * Alternatively, provide a file path to a `.gmt` or `.gmx` file to import
#'     custom gene sets.
#' - `data.frame`: A data frame with 1–4 columns. Expected columns include:
#'   * (Optional) Gene set ID
#'   * (Optional) term name
#'   * (Optional) description
#'   * Gene symbols or features
#' - `list`: A named or unnamed list where each element is a character vector of
#'   genes. Optionally, attach `term` or `description` attributes to each
#'   element, or `terms` / `descriptions` attributes to the entire list.
#' - `OrgDb`: An
#'   [`OrgDb`](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb)
#'   object (e.g., `org.Hs.eg.db`) to extract GO gene sets.
#' @param ... Additional arguments passed to the method.
#' @inheritParams gs_filter
#' @param _arg A single string of the argument name for `gs`, used by the
#' internal to provide better message.
#' @inheritParams keggdb
#' @return A `xbio_genesets` object.
#' @examples
#' # Read from a KEGG source
#' genesets("kegg")
#' repr_genesets("kegg", min_size = 5, max_size = 500)
#'
#' # Load from a local GMT file
#' \dontrun{
#' genesets("path/to/genesets.gmt")
#' repr_genesets("path/to/genesets.gmt", min_size = 5, max_size = 500)
#' }
#'
#' # Construct from a list
#' genesets(list(set1 = c("A", "B", "C"), set2 = c("X", "Y")))
#'
#' # From a data.frame
#' df <- data.frame(set = c("A", "A", "B"), gene = c("x", "y", "z"))
#' genesets(df)
#'
#' @export
repr_genesets <- function(gs, ..., min_size = NULL, max_size = NULL) {
    out <- genesets(gs, ...)
    out <- gs_lapply(gs, gs_clean.xbio_geneset)
    if (!all(keep <- list_sizes(out) > 0L)) {
        cli::cli_warn(paste(
            "Removing {sum(!keep)} invalid gene set{?s}",
            "(all are empty string or missing value)"
        ))
        out <- out[keep]
    }
    gs_filter(out, min_size = min_size, max_size = max_size)
}

#' @export
#' @rdname repr_genesets
genesets <- function(gs, ..., `_arg` = NULL) {
    genesets <- function(gs, ..., `_arg` = NULL) {
        if (missing(gs)) {
            cli::cli_abort("{.arg {`_arg`}} must be provided")
        }
        UseMethod("genesets")
    }
    # Tricks to do common work in S3 generic
    # Define the default argument of `_arg`
    genesets(gs, ..., `_arg` = `_arg` %||% "gs")
}

#' @export
genesets.default <- function(gs, ..., `_arg` = NULL) {
    cli::cli_abort("Cannot extract genesets from {.obj_type_friendly {gs}}")
}

#' @export
genesets.xbio_genesets <- function(gs, ..., `_arg` = NULL) gs

#' @param select A character or integer vector select the columns to be used for
#' a data.frame input.
#' @export
#' @rdname repr_genesets
genesets.data.frame <- function(gs, ..., select = NULL, `_arg` = NULL) {
    rlang::check_dots_empty()
    if (inherits(gs, "data.table")) gs <- as.data.frame(gs)
    if (!is.null(select)) {
        select <- vec_as_location(
            select,
            n = ncol(gs),
            names = names(gs),
            missing = "error"
        )
        if (length(select) == 0L || length(select) > 4L) {
            cli::cli_abort("{.arg select} must be of length [1, 4].")
        }
        gs <- gs[select]
    }
    vec_cast(gs, new_genesets())
}

#' @param ids A character vector of gene set IDs for list input. If `NULL`, the
#' names of the input list will be used.
#' @param terms A character vector of gene set terms for list input.
#' @param descriptions A character vector of gene set descriptions for a list
#' input.
#' @export
#' @rdname repr_genesets
genesets.list <- function(gs, ..., ids = NULL, terms = NULL,
                          descriptions = NULL, `_arg` = NULL) {
    rlang::check_dots_empty()
    if (!is.null(ids)) {
        ids <- vec_cast(ids, character())
        if (vec_size(gs) != vec_size(ids)) {
            cli::cli_abort(paste(
                "{.arg ids} must be",
                "the same length of the input list {.arg {`_arg`}}."
            ))
        }
        attr(gs, "ids") <- terms
    }
    if (!is.null(terms)) {
        terms <- vec_cast(terms, character())
        if (vec_size(gs) != vec_size(terms)) {
            cli::cli_abort(paste(
                "{.arg terms} must be",
                "the same length of the input list {.arg {`_arg`}}."
            ))
        }
        attr(gs, "terms") <- terms
    }
    if (!is.null(descriptions)) {
        descriptions <- vec_cast(descriptions, character())
        if (vec_size(gs) != vec_size(descriptions)) {
            cli::cli_abort(paste(
                "{.arg descriptions} must be",
                "the same length of the input list {.arg {`_arg`}}."
            ))
        }
        attr(gs, "descriptions") <- descriptions
    }
    vec_cast(gs, new_genesets())
}

#' @export
#' @rdname repr_genesets
genesets.character <- function(gs, ..., verbose = TRUE, `_arg` = NULL) {
    if (!.rlang_check_string(gs, allow_empty = FALSE)) {
        cli::cli_abort("{.arg {`_arg`}} must be a single non-empty string")
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
        orgdb <- rlang::try_fetch(
            getExportedValue(gs, gs),
            error = function(cnd) {
                cli::cli_abort(
                    "{.arg {`_arg`}} is not a valid {.cls OrgDb} package"
                )
            }
        )
        if (inherits(orgdb, "OrgDb")) {
            cli::cli_abort(
                "{.arg {`_arg`}} is not a valid {.cls OrgDb} package"
            )
        }
        return(genesets(gs = orgdb, ...))
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
#' @rdname repr_genesets
genesets.OrgDb <- function(gs, ..., keytype = "SYMBOL", ontology = NULL,
                           `_arg` = NULL) {
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
    out <- vec_slice(out, !is.na(out$GOALL))
    if (!is.null(ontology)) {
        out <- vec_slice(out, out$ONTOLOGYALL %in% ontology)
    }
    genesets(out, select = c("GOALL", "ONTOLOGYALL", keytype))
}

#' @export
genesets.xbio_kegg_genesets <- function(gs, ..., `_arg` = NULL) {
    vec_cast(gs, new_genesets(), x_arg = `_arg`)
}

genesets_kegg <- function(link = NULL, database = NULL, strategy = NULL,
                          save = NULL, verbose = TRUE, ...) {
    rlang::check_dots_empty()
    kegg_genesets(
        link = link,
        database = database,
        strategy = strategy,
        save = save,
        verbose = verbose
    )
}

genesets_go <- function(organism = NULL, strategy = NULL, save = NULL,
                        verbose = TRUE, ...) {
    rlang::check_dots_empty()
    godb(
        database = "genesets", organism = organism,
        strategy = strategy, save = save, verbose = verbose
    )
}
