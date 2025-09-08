#' Extract attributes from a geneset or genesets object
#'
#' These functions extract specific attributes from [geneset()] or
#' [genesets()] objects:
#'  - IDs (`gs_ids()`)
#'  - terms (`gs_terms()`)
#'  - descriptions (`gs_descs()`).
#'
#' @param gs A [geneset()] or [genesets()] object.
#'
#' @return A character vector containing the extracted attributes.
#' If `gs` is a `geneset` object, a character scalar is returned.
#'
#' @export
gs_ids <- function(gs) UseMethod("gs_ids")

#' @export
gs_ids.xbio_geneset <- function(gs) {
    o <- attr(gs, "id", exact = TRUE)
    if (is.null(o) || !is.character(o) || length(o) != 1L) {
        NA_character_
    } else {
        o
    }
}

#' @export
gs_ids.xbio_genesets <- function(gs) {
    gs_vapply(gs, gs_ids.xbio_geneset, character(1L))
}

#' @export
#' @rdname gs_ids
gs_terms <- function(gs) UseMethod("gs_terms")

#' @export
gs_terms.xbio_geneset <- function(gs) {
    o <- attr(gs, "term", exact = TRUE)
    if (is.null(o) || !is.character(o) || length(o) != 1L) {
        NA_character_
    } else {
        o
    }
}

#' @export
gs_terms.xbio_genesets <- function(gs) {
    gs_vapply(gs, gs_terms.xbio_geneset, character(1L))
}

#' @export
#' @rdname gs_ids
gs_descs <- function(gs) UseMethod("gs_descs")

#' @export
gs_descs.xbio_geneset <- function(gs) {
    o <- attr(gs, "description", exact = TRUE)
    if (is.null(o) || !is.character(o) || length(o) != 1L) {
        NA_character_
    } else {
        o
    }
}

#' @export
gs_descs.xbio_genesets <- function(gs) {
    gs_vapply(gs, gs_descs.xbio_geneset, character(1L))
}

#######################################################
#' Map a geneset or genesets object to a different keytype
#'
#' This function maps the gene identifiers in a [`geneset()`] or [`genesets()`]
#' object from one keytype to another, using either an
#' [`AnnotationDb`][AnnotationDbi::mapIds] or a [`Mart`][biomaRt::useMart]
#' object.
#'
#' @inheritParams gs_ids
#' @param annodb An `AnnotationDb` object or a character string naming the
#' annotation package.
#' @param mart A [`Mart`][biomaRt::useMart] object (for use with
#' `gs_biomart()`).
#' @param key_source A string indicating the keytype used in `gs`.
#' @param key_target A string indicating the desired target keytype.
#' @inheritParams AnnotationDbi::mapIds
#' @param ... Additional arguments passed to
#' [`mapIds()`][AnnotationDbi::mapIds] or [`getBM()`][biomaRt::getBM].
#'
#' @return An object of the same class as `gs`, with mapped gene identifiers.
#' @export
gs_map <- function(gs, annodb, key_source, key_target,
                   multiVals = "first", ...) {
    assert_string(key_source, allow_empty = FALSE)
    assert_string(key_target, allow_empty = FALSE)
    # Infer the annodb ----------------------------
    if (is.character(annodb)) {
        check_bioc_installed(annodb)
        annodb <- rlang::try_fetch(
            getExportedValue(annodb, annodb),
            error = function(cnd) {
                cli::cli_abort(
                    "{.arg annodb} is not a valid {.cls AnnotationDb} package"
                )
            }
        )
    }
    if (!inherits(annodb, "AnnotationDb")) {
        cli::cli_abort("{.arg annodb} must be a {.cls AnnotationDb} object")
    }
    gs_map <- function(gs, annodb, key_source, key_target, multiVals, ...) {
        UseMethod("gs_map")
    }
    gs_map(gs, annodb, key_source, key_target, multiVals, ...)
}

#' @export
gs_map.xbio_geneset <- function(gs, annodb, key_source, key_target,
                                multiVals = "first", ...) {
    if (vec_size(gs) == 0L) return(gs) # styler: off
    out <- AnnotationDbi::mapIds(
        x = annodb,
        keys = as.character(gs),
        column = key_target,
        keytype = key_source,
        multiVals = multiVals,
        ...
    )
    out <- as.character(out) # out can be a list
    out[is.na(out) | out == ""] <- NA_character_
    vec_restore(out, gs)
}

#' @export
gs_map.xbio_genesets <- function(gs, annodb, key_source, key_target,
                                 multiVals = "first", ...) {
    if (vec_size(gs) == 0L) return(gs) # styler: off
    # we'll filter multiVals manually
    if (rlang::is_string(multiVals, "filter")) {
        filter <- TRUE
        multiVals <- "list"
    } else {
        filter <- FALSE
    }
    # mapping the the genes in genesets into key_target
    mapped <- AnnotationDbi::mapIds(
        x = annodb,
        keys = list_unchop(lapply(vec_data(gs), as.character)),
        column = key_target,
        keytype = key_source,
        multiVals = multiVals,
        ...
    )
    mapped <- as.list(mapped)
    out <- .mapply(function(new, old) {
        if (filter) new <- new[list_sizes(new) == 1L]
        new <- as.character(new)
        new[is.na(new) | new == ""] <- NA_character_
        vec_restore(new, old)
    }, list(new = vec_chop(mapped, sizes = list_sizes(gs)), old = gs), NULL)
    vec_restore(out, gs)
}

#' @param mart A [`Mart`][biomaRt::useMart] object.
#' @export
#' @rdname gs_map
gs_biomart <- function(gs, mart, key_source, key_target, ...) {
    assert_string(key_source, allow_empty = FALSE)
    assert_string(key_target, allow_empty = FALSE)
    check_bioc_installed("biomaRt")
    if (missing(mart) || !inherits(mart, "Mart")) {
        cli::cli_abort(c(
            "{.arg mart} must be a valid {.cls Mart} object.",
            i = "To create a {.cls Mart} object use the function: {.fn biomaRt::useMart}",
            i = "Check {.code ?biomaRt::useMart} for more information."
        ))
    }
    UseMethod("gs_biomart")
}

#' @export
gs_biomart.xbio_geneset <- function(gs, mart, key_source, key_target, ...) {
    if (vec_size(gs) == 0L) return(gs) # styler: off
    out <- biomaRt::getBM(
        mart = mart,
        values = as.character(gs),
        attributes = key_target,
        filters = key_source,
        ...
    )
    out <- as.character(.subset2(out, key_target))
    out[is.na(out) | out == ""] <- NA_character_
    vec_restore(out, gs)
}

#' @export
gs_biomart.xbio_genesets <- function(gs, ...) {
    if (vec_size(gs) == 0L) return(gs) # styler: off
    # mapping the the genes in genesets into keytype
    gs_lapply(gs, gs_biomart.xbio_geneset, ...)
}

###########################################################
#' Filter genesets by size
#'
#' This function filters a [`genesets()`] object by removing gene sets that do
#' not fall within the specified size range.
#'
#' @param gs A [`genesets()`] object (of class `"xbio_genesets"`).
#' @param min_size An integer specifying the minimum number of genes
#' a gene set must contain to be retained. If `NULL`, no lower bound is applied.
#' @param max_size An integer specifying the maximum number of genes
#' a gene set may contain to be retained. If `NULL`, no upper bound is applied.
#'
#' @details
#' Gene sets whose sizes fall outside the interval defined by
#' `min_size` and `max_size` will be removed. A message is printed
#' indicating how many sets were removed, if any.
#'
#' @return A filtered [`genesets()`] object containing only gene sets
#' within the specified size range.
#'
#' @export
gs_filter <- function(gs, min_size = NULL, max_size = NULL) {
    assert_s3_class(gs, "xbio_genesets")
    assert_number_whole(min_size, min = 0, allow_null = TRUE)
    assert_number_whole(max_size, min = 0, allow_null = TRUE)
    if (is.null(min_size) && is.null(max_size)) {
        return(gs)
    }
    sizes <- list_sizes(gs)
    if (!is.null(min_size) && !is.null(max_size)) {
        keep <- sizes >= min_size & sizes <= max_size
        out_pattern <- sprintf("[%d, %d]", min_size, max_size)
    } else if (!is.null(min_size)) {
        keep <- sizes >= min_size
        out_pattern <- sprintf("[%d, Inf)", min_size)
    } else {
        keep <- sizes <= max_size
        out_pattern <- sprintf("(0, %d]", max_size)
    }
    if (!all(keep)) {
        cli::cli_inform(c(
            ">" = sprintf(
                "Removing {sum(!keep)} gene set{?s} with size out of %s",
                out_pattern
            )
        ))
        gs <- gs[keep]
    }
    gs
}

#' Clean gene sets by removing empty or invalid entries
#'
#' `gs_clean()` removes missing values and empty strings from each gene set.
#' Gene sets that become empty after cleaning are discarded.
#'
#' @inheritParams gs_ids
#' @return An object of the same class as `gs`, with invalid entries and empty
#' gene sets removed.
#' @examples
#' gs <- genesets(list(
#'     set1 = c("A", "B", NA, ""),
#'     set2 = c("", NA),
#'     set3 = c("C", "C")
#' ))
#' gs_clean(gs)
#'
#' @export
gs_clean <- function(gs) UseMethod("gs_clean")

#' @export
gs_clean.xbio_geneset <- function(gs) {
    gs <- vec_unique(gs)
    gs[!is.na(gs) & as.character(gs) != ""]
}

#' @export
gs_clean.xbio_genesets <- function(gs, ...) {
    out <- gs_lapply(gs, gs_clean.xbio_geneset)
    if (!all(keep <- list_sizes(out) > 0L)) {
        cli::cli_inform(c(
            ">" = paste(
                "Removing {sum(!keep)} invalid gene set{?s}",
                "(all are empty string or missing value)"
            )
        ), class = "gs_clean_message")
        out <- out[keep]
    }
    out
}

gs_lapply <- function(gs, ...) vec_restore(lapply(vec_data(gs), ...), gs)

gs_vapply <- function(gs, ...) vapply(vec_data(gs), ..., USE.NAMES = FALSE)
