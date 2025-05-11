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

new_geneset <- function(id = NA_character_, geneset = character(),
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

#' @export
vec_ptype2.GeneSet.NULL <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.NULL.GeneSet <- function(x, y, ...) {
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

#' @export
vec_cast.xbio_geneset.GeneSet <- function(x, to, ...,
                                          x_arg = caller_arg(x),
                                          to_arg = "",
                                          call = caller_env()) {
    slots <- methods::slotNames(x)
    values <- lapply(slots, function(nm) methods::slot(x, nm))
    names(values) <- slots
    rlang::inject(geneset(
        values$setName, values$geneIds,
        term = values$shortDescription,
        description = values$longDescription,
        !!!values[
            !names(values) %in% c(
                "setName", "geneIds",
                "shortDescription", "longDescription"
            )
        ]
    ))
}

#' @export
vec_cast.GeneSet.xbio_geneset <- function(x, to, ...,
                                          x_arg = caller_arg(x),
                                          to_arg = "",
                                          call = caller_env()) {
    slots <- c(
        "geneIdType", "setIdentifier", "organism",
        "pubMedIds", "urls", "contributor", "version", "creationDate",
        "collectionType"
    )
    names(slots) <- slots
    slots <- lapply(slots, function(slot) attr(x, slot, exact = TRUE))
    slots <- slots[vapply(slots, is.null, logical(1L), USE.NAMES = FALSE)]
    rlang::inject(GSEABase::GeneSet(
        as.character(x),
        setName = attr(x, "id"),
        shortDescription = attr(x, "term"),
        longDescription = attr(x, "description"),
        !!!slots
    ))
}
