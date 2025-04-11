genesets <- function(x, ...) UseMethod("genesets")

filter_genesets <- function(x, min_size, max_size) {
    if (is.null(min_size) && is.null(max_size)) {
        return(x)
    }
    sizes <- lengths(x)
    if (!is.null(min_size) && !is.null(max_size)) {
        x[sizes >= min_size & sizes <= max_size]
    } else if (!is.null(min_size)) {
        x[sizes >= min_size]
    } else {
        x[sizes <= max_size]
    }
}

#' @export
genesets.character <- function(x, ...) {

}
