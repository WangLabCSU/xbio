genesets <- function(x, ...) UseMethod("genesets")

filter_genesets <- function(x, min_size, max_size) {
    sizes <- lengths(x)
    x[sizes >= min_size & sizes <= max_size]
}

#' @export
genesets.character <- function(x, ...) {

}
