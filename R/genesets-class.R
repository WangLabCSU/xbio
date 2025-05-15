#' @importFrom rlang caller_arg
new_genesets <- function(genesets = list(),
                         ids = NULL, terms = NULL, descriptions = NULL, ...,
                         arg_genesets = caller_arg(genesets),
                         arg_ids = caller_arg(ids),
                         arg_terms = caller_arg(terms),
                         arg_descriptions = caller_arg(descriptions)) {
    if (is.null(ids)) {
        ids <- names(genesets)
    } else {
        ids <- vec_cast(ids, character(), x_arg = arg_ids)
        if (vec_size(ids) != vec_size(genesets)) {
            cli::cli_abort(paste(
                "{.arg {arg_ids}} ({vec_size(ids)}) must have",
                "the same length of {.arg {arg_genesets}} ({vec_size(genesets)})"
            ))
        }
    }
    if (is.null(terms)) {
        terms <- rep_len(list(NULL), vec_size(genesets))
    } else {
        terms <- vec_cast(terms, character(), x_arg = arg_terms)
        if (vec_size(terms) != vec_size(genesets)) {
            cli::cli_abort(paste(
                "{.arg {arg_terms}} ({vec_size(terms)}) must have",
                "the same length of {.arg {arg_genesets}} ({vec_size(genesets)})"
            ))
        }
    }
    if (is.null(descriptions)) {
        descriptions <- rep_len(list(NULL), vec_size(genesets))
    } else {
        descriptions <- vec_cast(
            descriptions, character(),
            x_arg = arg_descriptions
        )
        if (vec_size(descriptions) != vec_size(genesets)) {
            cli::cli_abort(paste(
                "{.arg {arg_descriptions}} ({vec_size(descriptions)}) must have",
                "the same length of {.arg {arg_genesets}} ({vec_size(genesets)})"
            ))
        }
    }
    genesets <- .mapply(
        new_geneset, list(
            id = ids %||% rep_len(list(NULL), vec_size(genesets)),
            geneset = genesets, term = terms, description = descriptions
        ),
        NULL
    )
    names(genesets) <- ids
    new_vctr(genesets, ..., class = "xbio_genesets")
}

S3_genesets <- new_S3_class("xbio_genesets") # Used by S7

#' @export
`names<-.xbio_genesets` <- function(x, value) {
    if (!is.null(value)) {
        value <- vec_cast(value, character(), x_arg = "names")
        vec_check_size(value, vec_size(x), arg = "names")
    }
    # we always ensure the names is the same with the element id
    old <- x
    x <- .mapply(
        function(geneset, id) {
            attr(geneset, "id") <- id
            geneset
        },
        list(geneset = x, id = value %||% rep_len(list(NULL), vec_size(x))),
        NULL
    )
    vec_restore(NextMethod(), old)
}

#' @export
vec_ptype_abbr.xbio_genesets <- function(x, ...) "genesets"

#' @export
obj_print_header.xbio_genesets <- function(x, ...) {
    cli::cat_line("<", "genesets", "[", vec_size(x), "]>")
    invisible(x)
}

#' @export
obj_print_data.xbio_genesets <- function(x, geneset_trunc = 6L,
                                         trunc = 6L, ...) {
    size <- vec_size(x)
    if (size == 0L) return(invisible(x)) # styler: off
    trunc <- max(as.integer(trunc), 2L)
    ids <- gs_ids(x)
    if (size > trunc) {
        before <- ceiling(trunc / 2L)
        after <- rev(seq_len(trunc - before))
        before <- seq_len(before)
        x <- vec_c(vec_slice(x, before), vec_slice(x, size - after))
        ids <- vec_c(vec_slice(ids, before), vec_slice(ids, size - after))
    }
    geneset_trunc <- max(as.integer(geneset_trunc), 2L)
    geneset_before <- ceiling(geneset_trunc / 2L)
    geneset_after <- rev(seq_len(geneset_trunc - geneset_before))
    geneset_before <- seq_len(geneset_before)
    content <- vapply(x, function(geneset) {
        geneset <- vec_cast(geneset, character())
        if (vec_size(geneset) > geneset_trunc) {
            geneset <- vec_c(
                vec_slice(geneset, geneset_before), "...",
                vec_slice(geneset, vec_size(geneset) - geneset_after)
            )
        }
        paste0(geneset, collapse = ", ")
    }, character(1L), USE.NAMES = FALSE)
    output <- paste(
        format(ids, justify = "right"),
        format(content, justify = "left"),
        sep = ": "
    )
    if (size > trunc) output <- vec_c(output[before], "...", output[-before])
    cat(output, sep = "\n")
}

##################################################
#' @export
vec_ptype2.xbio_genesets.xbio_genesets <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.xbio_genesets.NULL <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.NULL.xbio_genesets <- function(x, y, ...) {
    y
}

#' @export
vec_ptype2.xbio_genesets.data.frame <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.data.frame.xbio_genesets <- function(x, y, ...) {
    y
}

#' @export
vec_ptype2.xbio_genesets.list <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.list.xbio_genesets <- function(x, y, ...) {
    y
}

##################################################
#' @export
vec_cast.xbio_genesets.xbio_genesets <- function(x, to, ...) {
    x
}

#' @export
vec_cast.data.frame.xbio_genesets <- function(x, to, ...) {
    ids <- gs_ids(x)
    terms <- gs_terms(x)
    descriptions <- gs_descs(x)
    genesets <- lapply(x, vec_cast, to = character())
    new_data_frame(list(
        ids = ids,
        terms = terms,
        descriptions = descriptions,
        genesets = genesets
    ))
}

#' @export
vec_cast.xbio_genesets.data.frame <- function(x, to, ...) {
    # Nest the term and descriptions columns
    if (ncol(x) == 1L) {
        new_genesets(list(x[[1L]]))
    } else if (ncol(x) == 2L) {
        ids <- vec_cast(x[[1L]], character(), x_arg = "the 1st column")
        locs <- vec_group_loc(ids)
        new_genesets(
            vec_chop(x[[2L]], .subset2(locs, "loc")),
            ids = .subset2(locs, "key")
        )
    } else if (ncol(x) == 3L) {
        ids <- vec_cast(x[[1L]], character(), x_arg = "the 1st column")
        descriptions <- vec_cast(x[[2L]], character(), x_arg = "the 2nd column")
        locs <- vec_group_loc(data_frame(
            ids = ids,
            descriptions = descriptions
        ))
        keys <- .subset2(locs, "key")
        if (all(is.na(descriptions <- .subset2(keys, "descriptions")))) {
            descriptions <- NULL
        }
        new_genesets(
            vec_chop(x[[3L]], .subset2(locs, "loc")),
            ids = .subset2(keys, "ids"),
            descriptions = descriptions
        )
    } else if (ncol(x) == 4L) {
        ids <- vec_cast(x[[1L]], character(), x_arg = "the 1st column")
        terms <- vec_cast(x[[2L]], character(), x_arg = "the 2nd column")
        descriptions <- vec_cast(x[[3L]], character(), x_arg = "the 3nd column")
        locs <- vec_group_loc(data_frame(
            ids = ids, terms = terms,
            descriptions = descriptions
        ))
        keys <- .subset2(locs, "key")
        if (all(is.na(terms <- .subset2(keys, "terms")))) {
            terms <- NULL
        }
        if (all(is.na(descriptions <- .subset2(keys, "descriptions")))) {
            descriptions <- NULL
        }
        new_genesets(
            vec_chop(x[[4L]], .subset2(locs, "loc")),
            ids = .subset2(keys, "ids"),
            terms = terms, descriptions = descriptions
        )
    } else {
        cli::cli_abort(paste(
            "Conversion to {.cls genesets} require a data frame",
            "of 2, 3, or 4 columns"
        ))
    }
}

#' @export
vec_cast.xbio_genesets.xbio_kegg_genesets <- function(x, to, ...) {
    vec_cast.xbio_genesets.data.frame(
        new_data_frame(
            list(
                ids = .subset2(x, "ids"),
                terms = .subset2(x, "terms"),
                descriptions = rep_len(NA_character_, vec_size(x)),
                genesets = .subset2(x, "features")
            )
        ),
        to, ...
    )
}

#' @export
as.data.frame.xbio_genesets <- function(x, ...) vec_cast(x, new_data_frame())

#' @export
vec_cast.list.xbio_genesets <- function(x, to, ...) vec_data(x)

#' @export
vec_cast.xbio_genesets.list <- function(x, to, ...) {
    if (is.null(terms <- attr(x, "terms"))) {
        terms <- vapply(x, function(e) {
            if (is.null(o <- attr(x, "term"))) {
                NA_character_
            } else {
                as.character(o)
            }
        }, character(1L), USE.NAMES = FALSE)
        if (all(is.na(terms))) terms <- NULL
    } else {
        terms <- as.character(terms)
        if (vec_size(x) != vec_size(terms)) {
            cli::cli_abort(paste(
                "attribute {.field terms} must be",
                "the same length of the input list."
            ))
        }
    }
    if (is.null(descriptions <- attr(x, "descriptions"))) {
        descriptions <- vapply(x, function(e) {
            if (is.null(o <- attr(x, "description"))) {
                NA_character_
            } else {
                as.character(o)
            }
        }, character(1L), USE.NAMES = FALSE)
        if (all(is.na(descriptions))) descriptions <- NULL
    } else {
        descriptions <- as.character(descriptions)
        if (vec_size(x) != vec_size(descriptions)) {
            cli::cli_abort(paste(
                "attribute {.field descriptions} must be",
                "the same length of the input list."
            ))
        }
    }
    new_genesets(x, terms = terms, descriptions = descriptions)
}

#' @export
as.list.xbio_genesets <- function(x, ...) vec_cast(x, list())

#' @export
vec_math.xbio_genesets <- function(.fn, .x, ...) {
    stop_incompatible_op("vec_math", x = .x)
}

#' @export
`length<-.xbio_genesets` <- function(x, value) {
    size <- vec_size(x)
    if (value > size) {
        cli::cli_abort(
            "Cannot set length greater than the current number of gene sets."
        )
    }
    vec_slice(x, seq_len(value))
}

#' @export
rep.xbio_genesets <- function(x, ...) {
    out <- NextMethod()
    vec_restore(out, x)
}
