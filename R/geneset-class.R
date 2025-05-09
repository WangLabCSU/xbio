#' Create a single geneset
#'
#' @param id A single string of the geneset ID.
#' @param geneset A character of geneset elements.
#' @param term A single string of the geneset term.
#' @param description A single string of the geneset description.
#' @param ... Additional attributes to be added to the geneset.
#' @export
geneset <- function(id, geneset, term = NULL, description = NULL, ...) {
    assert_string(id, allow_empty = FALSE)
    assert_string(term,
        allow_empty = TRUE,
        allow_na = TRUE,
        allow_null = TRUE
    )
    assert_string(description,
        allow_empty = TRUE,
        allow_na = TRUE,
        allow_null = TRUE
    )
    UseMethod("geneset")
}

#' @export
#' @rdname geneset
geneset.character <- function(id, geneset, term = NULL,
                              description = NULL, ...) {
    new_geneset(
        id = id, geneset = geneset,
        term = term, description = description, ...
    )
}

#' @export
#' @rdname geneset
geneset.default <- function(id, geneset, term = NULL,
                            description = NULL, ...) {
    geneset <- vec_cast(geneset, character())
    new_geneset(
        id = id, geneset = geneset,
        term = term, description = description, ...
    )
}

new_geneset <- function(id, geneset = character(),
                        term = NULL,
                        description = NULL, ...,
                        class = NULL) {
    geneset <- vec_unique(geneset)
    new_vctr(geneset,
        id = id,
        term = term,
        description = description,
        ...,
        class = c(class, "xbio_geneset")
    )
}

#' @export
obj_print_header.xbio_geneset <- function(x, ...) {
    cli::cat_line(
        "<(",
        attr(x, "id", exact = TRUE),
        ") geneset",
        "[", vec_size(x), "]>"
    )
    invisible(x)
}

##################################################
#' @export
vec_ptype2.xbio_geneset.NULL <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.NULL.xbio_geneset <- function(x, y, ...) {
    y
}

# WE ALWAYS RETRUN A CHARACTER WHEN COMBINED WITH CHARACTER
# Since the elements of a geneset must be unique, if we return
# a geneset, it's not possible to keep only unique elements
#' @export
vec_ptype2.xbio_geneset.character <- function(x, y, ...) {
    y
}

#' @export
vec_ptype2.character.xbio_geneset <- function(x, y, ...) {
    x
}

#############################################################
#' @export
vec_cast.character.xbio_geneset <- function(x, to, ...,
                                            x_arg = caller_arg(x),
                                            to_arg = "",
                                            call = caller_env()) {
    vec_data(x)
}
