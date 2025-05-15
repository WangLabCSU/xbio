#' Extract matrix
#'
#' @param matrix Input object used to extract the matrix.
#' @inheritParams repr_genesets
#' @return A matrix object.
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

#' @export
repr_matrix.default <- function(matrix, ..., `_arg` = NULL) {
    cli::cli_abort("Cannot extract matrix from {.obj_type_friendly {matrix}}")
}

#' @export
#' @rdname repr_matrix
repr_matrix.matrix <- function(matrix, ..., `_arg` = NULL) matrix

#' @export
#' @rdname repr_matrix
repr_matrix.dgCMatrix <- function(matrix, ..., `_arg` = NULL) matrix

#' @export
#' @rdname repr_matrix
repr_matrix.Seurat <- function(matrix, assay = NULL, ...,
                               dimred = NULL, n_dimred = NULL,
                               layer = "data", `_arg` = NULL) {
    if (!is.null(dimred)) {
        # value is expected to be a matrix or matrix-like object with number of
        # rows equal to ncol(x).
        mat <- getExportedValue("SeuratObject", "Embeddings")(
            matrix, reduction = dimred
        )
        if (!is.null(n_dimred)) {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
        t(mat)
    } else {
        getExportedValue("SeuratObject", "GetAssayData")(
            matrix, assay = assay, layer = layer
        )
    }
}

#' @export
#' @rdname repr_matrix
repr_matrix.SummarizedExperiment <- function(matrix, assay = NULL, ...,
                                             `_arg` = NULL) {
    getExportedValue("SummarizedExperiment", "assay")(matrix, assay %||% 1L)
}

#' @export
#' @rdname repr_matrix
repr_matrix.SingleCellExperiment <- function(matrix, assay = NULL, ...,
                                             dimred = NULL, n_dimred = NULL,
                                             `_arg` = NULL) {
    if (!is.null(dimred)) {
        # value is expected to be a matrix or matrix-like object with number of
        # rows equal to ncol(x).
        mat <- getExportedValue("SingleCellExperiment", "reducedDim")(
            matrix, dimred
        )
        if (!is.null(n_dimred)) {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
        t(mat)
    } else {
        getExportedValue("SummarizedExperiment", "assay")(matrix, assay %||% 1L)
    }
}
