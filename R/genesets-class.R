#' @importFrom rlang caller_arg
new_genesets <- function(genesets, ..., terms = NULL, descriptions = NULL,
                         arg_genesets = caller_arg(genesets),
                         arg_terms = caller_arg(terms),
                         arg_descriptions = caller_arg(descriptions)) {
    if (is.null(terms)) {
        if (!rlang::is_named2(genesets)) {
            cli::cli_abort(paste(
                "{.arg {arg_terms}} must be provided or {.arg {arg_genesets}}",
                "should be named"
            ))
        }
    } else {
        terms <- vec_cast(terms, character(), x_arg = arg_terms)
        if (vec_size(terms) != vec_size(genesets)) {
            cli::cli_abort(paste(
                "{.arg {arg_terms}} ({vec_size(terms)}) must have",
                "the same length of {.arg {arg_genesets}} ({vec_size(genesets)})"
            ))
        }
        names(genesets) <- terms
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
    }
    # we use new class to ensure `descriptions` are parallel with the data
    new_vctr(
        genesets,
        descriptions = descriptions,
        ...,
        class = "enricher_genesets"
    )
}

#' @export
`names<-.enricher_genesets` <- function(x, value) {
    if (is.null(value)) {
        cli::cli_abort("Cannot remove the names of genesets")
    }
    value <- vec_cast(value, character(), x_arg = "names")
    if (vec_any_missing(value) || any(value == "")) {
        cli::cli_abort("Names cannot be missing or empty.")
    }
    NextMethod()
}

#' @export
vec_proxy.enricher_genesets <- function(x, ...) {
    descriptions <- attr(x, "descriptions", exact = TRUE)
    terms <- names(x)
    genesets <- unclass(x)
    attributes(genesets) <- NULL
    if (is.null(descriptions)) {
        new_data_frame(list(terms = terms, genesets = genesets))
    } else {
        new_data_frame(
            list(
                terms = terms,
                descriptions = descriptions,
                genesets = genesets
            )
        )
    }
}

#' @export
vec_restore.enricher_genesets <- function(x, to, ...) {
    # restore the terms and descriptions from the proxy
    terms <- .subset2(x, "terms")
    attr(to, "descriptions") <- .subset2(x, "descriptions")

    # restore the attributes for the genesets
    x <- .subset2(x, "genesets")
    out <- NextMethod()
    names(out) <- terms
    out
}

#' @export
obj_print_data.enricher_genesets <- function(x, ...) {
    data <- vec_data(x)
    size <- vec_size(data)
    if (size > 6L) {
        data <- vec_c(vec_slice(data, 1:3), vec_slice(data, size - (3:1)))
    }
    content <- vapply(
        .subset2(data, "genesets"),
        function(geneset) {
            if (vec_size(geneset) > 6L) {
                geneset <- vec_c(
                    vec_slice(geneset, 1:3), "...",
                    vec_slice(geneset, vec_size(geneset) - (3:1))
                )
            }
            paste0(geneset, collapse = ", ")
        },
        character(1L),
        USE.NAMES = FALSE
    )
    output <- paste(
        format(.subset2(data, "terms"), justify = "right"),
        format(content, justify = "left"),
        sep = ": "
    )
    if (size > 6L) output <- vec_c(output[1:3], "...", output[4:6])
    cat(output, sep = "\n")
}

#' @export
vec_cast.enricher_genesets.enricher_genesets <- function(x, to, ...) {
    x
}

#' @export
vec_cast.data.frame.enricher_genesets <- function(x, to, ...) {
    vec_proxy(x)
}

#' @export
vec_cast.enricher_genesets.data.frame <- function(x, to, ...) {
    if (ncol(x) == 2L) {
        new_genesets(x[[2L]],
            terms = vec_cast(x[[1L]], character(), x_arg = "the 1st column")
        )
    } else if (ncol(x) == 3L) {
        new_genesets(x[[3L]],
            terms = vec_cast(x[[1L]], character(), x_arg = "the 1st column"),
            descriptions = vec_cast(
                x[[2L]], character(),
                x_arg = "the 2nd column"
            )
        )
    } else {
        cli::cli_abort("{.arg x} must be a data frame of 2 or 3 columns")
    }
}

#' @export
as.data.frame.enricher_genesets <- function(x, ...) {
    vec_cast(x, new_data_frame())
}

#' @export
vec_cast.list.enricher_genesets <- function(x, to, ...) unclass(x)

#' @export
vec_cast.enricher_genesets.list <- function(x, to, ...) {
    if (!rlang::is_named2(x)) {
        cli::cli_abort("{.arg x} must be a named list")
    }
    if (is.null(descriptions <- attr(x, "descriptions"))) {
        descriptions <- vapply(
            x,
            function(e) as.character(attr(x, "description")) %||% NA_character_,
            character(1L),
            USE.NAMES = FALSE
        )
        if (all(is.na(descriptions))) descriptions <- NULL
    } else {
        descriptions <- as.character(descriptions)
        if (vec_size(x) != vec_size(descriptions)) {
            descriptions <- NULL
        }
    }
    new_genesets(x, descriptions = descriptions)
}

#' @export
as.list.enricher_genesets <- function(x, ...) vec_cast(x, list())

#' @export
vec_math.enricher_genesets <- function(.fn, .x, ...) {
    stop_incompatible_op("vec_math", x = .x)
}

#' @export
`[.enricher_genesets` <- function(x, i, ...) {
    if (!missing(...)) {
        cli::cli_abort(
            "Can't index {.cls genesets} on dimensions greater than 1."
        )
    }
    data <- new_data_frame(vec_proxy(x), row.names = names(x))
    vec_restore(vec_slice(data, rlang::maybe_missing(i)), x)
}

#' @export
`[[.enricher_genesets` <- function(x, i, ...) {
    i <- vec_as_location2(i, n = length(x), names = names(x))
    out <- do.call(`[[`, list(vec_cast(x, list()), i, ...))
    descriptions <- attr(x, "descriptions", exact = TRUE)
    if (!is.null(descriptions)) {
        attr(out, "description") <- descriptions[i]
    }
    out
}

#' @export
`$.enricher_genesets` <- `[[.enricher_genesets`

#' @export
`length<-.enricher_genesets` <- function(x, value) {
    data <- vec_proxy(x)
    size <- vec_size(data)
    if (value > size) {
        cli::cli_abort(
            "Cannot set length greater than the current number of gene sets."
        )
    } else {
        i <- seq_len(value)
    }
    vec_restore(vec_slice(data, i), x)
}

#' @export
rep.enricher_genesets <- function(x, ...) {
    out <- lapply(vec_proxy(x), base::rep, ...)
    vec_restore(out, x)
}
