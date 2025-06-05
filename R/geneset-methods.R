#' Create a single geneset
#'
#' @param geneset A character of geneset elements.
#' @param ... Additional attributes to be added to the geneset.
#' @export
geneset <- function(geneset, ...) UseMethod("geneset")

#' @export
#' @rdname geneset
geneset.default <- function(geneset, id = NULL, term = NULL,
                            description = NULL, ...) {
    if (missing(geneset)) {
        geneset <- character()
    } else {
        geneset <- vec_cast(geneset, character())
    }
    geneset.character(
        geneset = geneset,
        id = id, term = term, description = description,
        ...
    )
}

#' @param id A single string of the geneset ID.
#' @param term A single string of the geneset term.
#' @param description A single string of the geneset description.
#' @export
#' @rdname geneset
geneset.character <- function(geneset, id = NULL, term = NULL,
                              description = NULL, ...) {
    assert_string(id,
        allow_empty = TRUE,
        allow_na = TRUE,
        allow_null = TRUE
    )
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
    new_geneset(
        geneset = geneset,
        id = id, term = term, description = description,
        ...
    )
}

#' @export
geneset.GeneSet <- function(geneset, ...) {
    rlang::check_dots_empty()
    vec_cast(geneset, new_geneset())
}
