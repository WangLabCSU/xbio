new_geneset <- function(geneset = character(),
                        term = NULL,
                        description = NULL, ...,
                        class = NULL) {
    geneset <- vec_unique(as.character(geneset))
    new_vctr(geneset,
        term = term,
        description = description,
        ...,
        class = c(class, "enricher_geneset")
    )
}

#' @export
obj_print_header.enricher_geneset <- function(x, ...) {
    cli::cat_line("<", "geneset", "[", vec_size(x), "]>")
    invisible(x)
}

##################################################
#' @export
vec_ptype2.enricher_geneset.NULL <- function(x, y, ...) {
    x
}

#' @export
vec_ptype2.NULL.enricher_geneset <- function(x, y, ...) {
    y
}

#############################################################
#' @export
vec_cast.character.enricher_geneset <- function(x, to, ...,
                                                x_arg = caller_arg(x),
                                                to_arg = "",
                                                call = caller_env()) {
    vec_data(x)
}
