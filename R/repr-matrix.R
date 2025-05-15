#' Extract matrix
#'
#' @param matrix Input object used to extract the matrix.
#' @inheritParams repr_genesets
#' @return A `xbio_matrix` object, which is a named numeric vector.
#' @export
repr_matrix <- function(matrix, ..., `_arg` = NULL) {
    repr_matrix <- function(matrix, ..., `_arg` = NULL) {
        if (missing(matrix)) {
            cli::cli_abort("{.arg {`_arg`}} must be provided")
        }
        UseMethod("repr_matrix")
    }
    # Tricks to do common work in S3 generic
    repr_matrix(matrix, ..., `_arg` = `_arg` %||% "matrix")
}

new_matrix <- function(x) structure(x, class = "xbio_matrix")

methods::setOldClass("xbio_matrix")

#' @export
repr_matrix.default <- function(matrix, ..., `_arg` = NULL) {
    cli::cli_abort("Cannot extract matrix from {.obj_type_friendly {matrix}}")
}

#' @export
repr_matrix.repr_matrix <- function(matrix, ..., `_arg` = NULL) matrix

#' @export
#' @rdname repr_matrix
repr_matrix.matrix <- function(matrix, ..., `_arg` = NULL) new_matrix(matrix)
