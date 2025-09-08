#' Create genesets representation
#'
#' @description
#' `repr_genesets()` is a wrapper around [`genesets()`] that additionally cleans
#' and filters the gene sets. It removes invalid entries, such as gene sets with
#' no genes or those containing only missing values. This function ensures the
#' returned object is clean and ready for downstream analysis.
#'
#' All analyses that depend on the gene set representation will use
#' `repr_genesets()` with its default arguments. If you want to customize
#' filtering, such as setting thresholds on gene set size (`min_size` /
#' `max_size`), you must call `repr_genesets()` manually.
#'
#' @inheritParams genesets
#' @inheritDotParams genesets
#' @return A `xbio_genesets` object.
#' @examples
#' # Read from a KEGG source
#' repr_genesets("kegg", min_size = 5, max_size = 500)
#'
#' # Load from a local GMT file
#' \dontrun{
#' repr_genesets("path/to/genesets.gmt", min_size = 5, max_size = 500)
#' }
#'
#' # Construct from a list
#' repr_genesets(
#'     list(
#'         set1 = c("A", "B", "C"),
#'         set2 = c("X", "Y")
#'     ),
#'     min_size = 3
#' )
#'
#' # From a data.frame
#' df <- data.frame(set = c("A", "A", "B"), gene = c("x", "y", "z"))
#' repr_genesets(df, min_size = 3)
#'
#' @export
repr_genesets <- function(gs, ..., min_size = NULL, max_size = NULL) {
    out <- genesets(gs, ...)
    out <- gs_clean(out)
    gs_filter(out, min_size = min_size, max_size = max_size)
}
