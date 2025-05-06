# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
write_gmt <- function(gs, path) {

}

read_gmt <- function(path, ...) {
    data <- strsplit(read_lines(path), "\t", fixed = TRUE)
    data <- data[lengths(data) > 2L]
    if (length(data) == 0L) cli::cli_abort("No genesets in {.path {path}}")
    ids <- vapply(data, .subset, 1L, character(1L), USE.NAMES = FALSE)
    descriptions <- vapply(data, .subset, 2L, character(1L), USE.NAMES = FALSE)
    genesets <- lapply(data, function(geneset) {
        trimws(.subset(geneset, -(1:2)), "both")
    })
    new_genesets(genesets, ids = ids, descriptions = descriptions)
}

# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29
read_gmx <- function(path, ...) {
    data <- utils::read.table(path, header = FALSE, sep = "\t")
    data <- data[lengths(data) > 2L]
    if (length(data) == 0L) cli::cli_abort("No genesets in {.path {path}}")
    ids <- vapply(data, .subset, 1L, character(1L), USE.NAMES = FALSE)
    descriptions <- vapply(data, .subset, 2L, character(1L), USE.NAMES = FALSE)
    genesets <- lapply(data, function(geneset) {
        trimws(.subset(geneset, -(1:2)), "both")
    })
    new_genesets(genesets, ids = ids, descriptions = descriptions)
}
