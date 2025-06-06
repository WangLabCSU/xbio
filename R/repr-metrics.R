#' Create ranking metrics representation
#'
#' @param metrics A named numeric ranked metrics.
#' @inheritParams repr_genesets
#' @return A named numeric vector with class `xbio_metrics`.
#' @export
repr_metrics <- function(metrics, ..., `_arg` = NULL) {
    repr_metrics <- function(metrics, ..., `_arg` = NULL) {
        if (missing(metrics)) {
            cli::cli_abort("{.arg {`_arg`}} must be provided")
        }
        UseMethod("repr_metrics")
    }
    # Tricks to do common work in S3 generic
    # Define the default argument of `_arg`
    repr_metrics(metrics, ..., `_arg` = `_arg` %||% "metrics")
}

new_metrics <- function(x) structure(x, class = "xbio_metrics")

is_valid_metrics <- function(metrics, names) {
    valid <- !vec_detect_missing(names) & names != ""
    if (!all(valid)) {
        n_wrong <- sum(!valid) # nolint
        cli::cli_warn("Removing {n_wrong} ranking metrics without valid names")
    }
    valid
}

#' @export
repr_metrics.default <- function(metrics, ..., `_arg` = NULL) {
    cli::cli_abort(
        "Cannot extract ranking metrics from {.obj_type_friendly {metrics}}"
    )
}

#' @export
repr_metrics.xbio_metrics <- function(metrics, ..., `_arg` = NULL) metrics

#########################################################
#' @export
#' @rdname repr_metrics
repr_metrics.numeric <- function(metrics, ..., `_arg` = NULL) {
    if (is.null(nms <- names(metrics))) {
        cli::cli_abort("{.arg {`_arg`}} must be a named numeric")
    }
    new_metrics(metrics[is_valid_metrics(metrics, nms)])
}

#' @export
#' @rdname repr_metrics
repr_metrics.data.frame <- function(metrics, ..., `_arg` = NULL) {
    if (ncol(metrics) < 2L) {
        cli::cli_abort(
            "{.arg {arg}} must be a data frame of at least 2 columns"
        )
    } else if (ncol(metrics) > 2L) {
        cli::cli_warn("Only the head 2 columns will be used")
    }
    nms <- vec_cast(metrics[[1L]], character(),
        x_arg = sprintf("the 1st column of %s", `_arg`)
    )
    metrics <- vec_cast(metrics[[2L]], numeric(),
        x_arg = sprintf("the 2nd column of %s", `_arg`)
    )
    names(metrics) <- nms
    new_metrics(metrics[is_valid_metrics(metrics, nms)])
}


#' @export
#' @rdname repr_metrics
repr_metrics.character <- function(metrics, ..., `_arg` = NULL) {
    if (!rlang::is_string(metrics) || !endsWith(tolower(metrics), ".rnk")) {
        cli::cli_abort(
            "{.arg {`_arg`}} must be a single string of {.field .rnk} file path"
        )
    }
    data <- read_rnk(metrics)

    # Call the data frame method
    repr_metrics(metrics = data, ..., `_arg` = `_arg`)
}


write_rnk <- function(prerank, path) {
    write_table(
        data.frame(names = names(prerank), value = prerank),
        path = path, col.names = FALSE, quote = FALSE
    )
}

read_rnk <- function(path) {
    read_table(path, header = FALSE, comment = "", quote = "")
}
