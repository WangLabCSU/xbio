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
    slots <- methods::slotNames(geneset)
    attrs <- lapply(slots, function(nm) methods::slot(geneset, nm))
    names(attrs) <- slots
    rlang::inject(geneset(
        attrs$setName, attrs$geneIds,
        term = attrs$shortDescription,
        description = attrs$longDescription,
        !!!attrs[
            !names(attrs) %in% c(
                "setName", "geneIds",
                "shortDescription", "longDescription"
            )
        ]
    ))
}

new_geneset <- function(geneset = character(),
                        id = NULL, term = NULL, description = NULL,
                        ...,
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
vec_ptype_abbr.xbio_geneset <- function(x, ...) "geneset"

#' @export
obj_print_header.xbio_geneset <- function(x, ...) {
    if (is.null(id <- attr(x, "id", exact = TRUE))) {
        cli::cat_line("<geneset", "[", vec_size(x), "]>")
    } else {
        cli::cat_line("<(", id, ") geneset[", vec_size(x), "]>")
    }
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

#' @export
vec_ptype2.xbio_geneset.character <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.character.xbio_geneset <- function(x, y, ...) {
    y
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
vec_cast.xbio_geneset.character <- function(x, to, ...,
                                            x_arg = caller_arg(x),
                                            to_arg = "",
                                            call = caller_env()) {
    geneset(x)
}

#' @export
vec_cast.xbio_geneset.GeneSet <- function(x, to, ...,
                                          x_arg = caller_arg(x),
                                          to_arg = "",
                                          call = caller_env()) {
    geneset(x)
}

#' @export
vec_cast.GeneSet.xbio_geneset <- function(x, to, ...,
                                          x_arg = caller_arg(x),
                                          to_arg = "",
                                          call = caller_env()) {
    check_bioc_installed("GSEABase")
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
        setName = attr(x, "id", exact = TRUE) %||% NA_character_,
        shortDescription = attr(x, "term", exact = TRUE) %||% NA_character_,
        longDescription = attr(x, "description", exact = TRUE) %||%
            NA_character_,
        !!!slots
    ))
}
