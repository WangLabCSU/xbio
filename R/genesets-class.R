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

validate_genesets <- function(gs) {
    if (vec_duplicate_any(names(gs))) {
        cli::cli_abort("the names of {.cls genesets} must be unique")
    }
    gs
}

#' @export
vec_proxy.enricher_genesets <- function(x, ...) {
    descriptions <- attr(x, "descriptions", exact = TRUE)
    terms <- names(x)
    genesets <- unclass(x)
    attributes(genesets) <- NULL
    if (is.null(descriptions)) {
        new_data_frame(
            list(terms = terms, genesets = genesets),
            row.names = terms
        )
    } else {
        new_data_frame(
            list(
                terms = terms,
                descriptions = descriptions,
                genesets = genesets
            ),
            row.names = terms
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
`[.enricher_genesets` <- function(x, i, ...) {
    if (!missing(...)) {
        cli::cli_abort(
            "Can't index {.cls genesets} on dimensions greater than 1."
        )
    }
    vec_restore(vec_slice(vec_proxy(x), rlang::maybe_missing(i)), x)
}

#' @export
`[[.enricher_genesets` <- function(x, i, ...) {
    nms <- names(x)
    x <- unclass(x)
    attributes(x) <- NULL
    names(x) <- nms
    x[[i]]
}

#' @export
`$.enricher_genesets` <- `[[.enricher_genesets`

#' @export
`length<-.enricher_genesets` <- function(x, value) {
    out <- vec_size_assign(vec_data(x), value)
    vec_restore(out, x)
}

#' @export
obj_print_data.enricher_genesets <- function(x, ...) {
    data <- vec_data(x)
    size <- vec_size(x)
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
rep.enricher_genesets <- function(x, ...) {
    out <- lapply(vec_data(x), base::rep, ...)
    vec_restore(out, x)
}
