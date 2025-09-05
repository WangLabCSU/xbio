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
as.character.xbio_geneset <- function(x) vec_cast(x, geneset())

#' @export
vec_cast.xbio_geneset.character <- function(x, to, ...,
                                            x_arg = caller_arg(x),
                                            to_arg = "",
                                            call = caller_env()) {
    new_geneset(x)
}

#' @export
vec_cast.xbio_geneset.GeneSet <- function(x, to, ...,
                                          x_arg = caller_arg(x),
                                          to_arg = "",
                                          call = caller_env()) {
    slots <- methods::slotNames(x)
    attrs <- lapply(slots, function(nm) methods::slot(x, nm))
    names(attrs) <- slots
    rlang::inject(geneset(
        attrs$setName,
        id = attrs$geneIds,
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
        setName = gs_ids(x),
        shortDescription = gs_terms(x),
        longDescription = gs_descs(x),
        !!!slots
    ))
}
