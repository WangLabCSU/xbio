new_genesets <- function(genesets, ..., terms = NULL, descriptions = NULL) {
    if (is.null(terms)) {
        if (!rlang::is_named(genesets)) {
            cli::cli_abort(paste(
                "{.arg terms} must be provided or {.arg genesets}",
                "should be named"
            ))
        }
        terms <- names(genesets)
    }
    new_vctr(
        genesets,
        terms = terms,
        descriptions = descriptions,
        ...,
        class = "enricher_genesets"
    )
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
    attr(to, "terms") <- .subset2(x, "terms")
    attr(to, "descriptions") <- .subset2(x, "descriptions")
    x <- .subset2(x, "genesets")
    NextMethod()
}

#' @export
names.enricher_genesets <- function(x) attr(x, "terms", exact = TRUE)

#' @export
`names<-.enricher_genesets` <- function(x, value) {
    if (is.null(value)) {
        cli::cli_abort("Cannot remove the names of genesets")
    } else {
        value <- as.character(value)
        vec_check_size(value, vec_size(x), arg = "names")
        attr(x, "terms") <- value
    }
    x
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

#' @method vec_cast enricher_genesets
#' @export
vec_cast.enricher_genesets <- function(x, to, ...) {
    UseMethod("vec_cast.enricher_genesets")
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
        new_genesets(.subset2(x, 2L), terms = .subset2(x, 1L))
    } else if (ncol(x) == 3L) {
        new_genesets(
            .subset2(x, 3L),
            terms = .subset2(x, 1L),
            descriptions = .subset2(x, 2L)
        )
    } else {
        cli::cli_abort("{.arg x} must be a data frame of 2 or 3 columns")
    }
}

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
    nms <- names(x)
    x <- unclass(x)
    names(x) <- nms
    out <- do.call("[[", list(x, i, ...))
    descriptions <- attr(x, "descriptions", exact = TRUE)
    if (!is.null(descriptions)) {
        names(descriptions) <- nms
        attr(out, "description") <- do.call("[[", list(descriptions, i))
    }
    out
}

#' @export
`$.enricher_genesets` <- function(x, i, ...) {
    nms <- names(x)
    x <- unclass(x)
    names(x) <- nms
    out <- do.call("$", list(x, i, ...))
    descriptions <- attr(x, "descriptions", exact = TRUE)
    if (!is.null(descriptions)) {
        names(descriptions) <- nms
        attr(out, "description") <- do.call("$", list(descriptions, i))
    }
    out
}

#' @export
`length<-.enricher_genesets` <- function(x, value) {
    out <- vec_size_assign(vec_proxy(x), value)
    vec_restore(out, x)
}

#' @export
rep.enricher_genesets <- function(x, ...) {
    out <- lapply(vec_proxy(x), base::rep, ...)
    vec_restore(out, x)
}
