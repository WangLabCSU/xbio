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
    if (all(is.na(ids))) ids <- NULL
    if (!is.null(terms)) {
        terms <- vec_cast(terms, character(), x_arg = arg_terms)
        if (vec_size(terms) != vec_size(genesets)) {
            cli::cli_abort(paste(
                "{.arg {arg_terms}} ({vec_size(terms)}) must have",
                "the same length of {.arg {arg_genesets}} ({vec_size(genesets)})"
            ))
        }
        if (all(is.na(terms))) terms <- NULL
    }
    if (!is.null(descriptions)) {
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
        if (all(is.na(descriptions))) descriptions <- NULL
    }
    genesets <- .mapply(
        new_geneset, list(
            geneset = genesets,
            id = ids %||% rep_len(list(NULL), vec_size(genesets)),
            term = terms %||% rep_len(list(NULL), vec_size(genesets)),
            description = descriptions %||%
                rep_len(list(NULL), vec_size(genesets))
        ),
        NULL
    )
    new_vctr(genesets, ..., class = "xbio_genesets")
}

S3_genesets <- new_S3_class("xbio_genesets") # Used by S7

#' @export
names.xbio_genesets <- function(x) gs_ids(x)

#' @export
`names<-.xbio_genesets` <- function(x, value) {
    if (!is.null(value)) {
        value <- vec_cast(value, character(), x_arg = "names")
        vec_check_size(value, vec_size(x), arg = "names")
    }
    # we always ensure the names is the same with the element id
    out <- .mapply(
        function(geneset, id) {
            attr(geneset, "id") <- id
            geneset
        },
        list(
            geneset = vec_data(x),
            id = value %||% rep_len(list(NULL), vec_size(x))
        ),
        NULL
    )
    vec_restore(out, x)
}

#' @export
`[.xbio_genesets` <- function(x, i, ...) {
    if (!missing(...)) {
        abort("Can't index genesets on dimensions greater than 1.")
    }
    index <- vec_as_location(i,
        n = vec_size(x), names = gs_ids(x),
        missing = "error"
    )
    vec_slice(x, index)
}

#' @export
`[[.xbio_genesets` <- function(x, i, ...) {
    index <- vec_as_location2(i,
        n = vec_size(x), names = gs_ids(x),
        missing = "error"
    )
    .subset2(x, index)
}

#' @export
`$.xbio_genesets` <- `[[.xbio_genesets`

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
    content <- gs_vapply(x, function(geneset) {
        geneset <- vec_cast(geneset, character())
        if (vec_size(geneset) > geneset_trunc) {
            geneset <- vec_c(
                vec_slice(geneset, geneset_before), "...",
                vec_slice(geneset, vec_size(geneset) - geneset_after)
            )
        }
        paste0(geneset, collapse = ", ")
    }, character(1L))
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
    genesets <- lapply(vec_data(x), vec_cast, to = character())
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
        new_genesets(
            vec_chop(x[[3L]], .subset2(locs, "loc")),
            ids = .subset2(keys, "ids"),
            descriptions = .subset2(keys, "descriptions")
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
        new_genesets(
            vec_chop(x[[4L]], .subset2(locs, "loc")),
            ids = .subset2(keys, "ids"),
            terms = .subset2(keys, "terms"),
            descriptions = .subset2(keys, "descriptions")
        )
    } else {
        cli::cli_abort(c(
            "Conversion to {.cls genesets} failed.",
            "x" = "A data frame with 1 to 4 columns is required."
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
vec_cast.list.xbio_genesets <- function(x, to, ...) {
    out <- vec_data(x)
    names(out) <- vapply(
        out, gs_ids.xbio_geneset,
        character(1L),
        USE.NAMES = FALSE
    )
    out
}

#' @export
vec_cast.xbio_genesets.list <- function(x, to, ...) {
    if (is.null(ids <- attr(x, "ids"))) {
        ids <- vapply(x, function(e) {
            if (is.null(o <- attr(x, "id")) || length(o) != 1L) {
                NA_character_
            } else {
                as.character(o)
            }
        }, character(1L), USE.NAMES = FALSE)
    } else {
        ids <- as.character(ids)
        if (vec_size(x) != vec_size(ids)) {
            cli::cli_warn(paste(
                "attribute {.field ids} must be",
                "the same length of the input list."
            ))
            ids <- NULL
        }
    }
    if (is.null(terms <- attr(x, "terms"))) {
        terms <- vapply(x, function(e) {
            if (is.null(o <- attr(x, "term")) || length(o) != 1L) {
                NA_character_
            } else {
                as.character(o)
            }
        }, character(1L), USE.NAMES = FALSE)
    } else {
        terms <- as.character(terms)
        if (vec_size(x) != vec_size(terms)) {
            cli::cli_warn(paste(
                "attribute {.field terms} must be",
                "the same length of the input list."
            ))
            terms <- NULL
        }
    }
    if (is.null(descriptions <- attr(x, "descriptions"))) {
        descriptions <- vapply(x, function(e) {
            if (is.null(o <- attr(x, "description")) || length(o) != 1L) {
                NA_character_
            } else {
                as.character(o)
            }
        }, character(1L), USE.NAMES = FALSE)
    } else {
        descriptions <- as.character(descriptions)
        if (vec_size(x) != vec_size(descriptions)) {
            cli::cli_warn(paste(
                "attribute {.field descriptions} must be",
                "the same length of the input list."
            ))
            descriptions <- NULL
        }
    }
    new_genesets(x, ids = ids, terms = terms, descriptions = descriptions)
}

#' @export
as.list.xbio_genesets <- function(x, ...) vec_cast(x, list())

#' @export
vec_math.xbio_genesets <- function(.fn, .x, ...) {
    stop_incompatible_op("vec_math", x = .x)
}

#' @export
`length<-.xbio_genesets` <- function(x, value) {
    if (value > vec_size(x)) {
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
