write_rnk <- function(prerank, path) {
    write_table(
        data.frame(names = names(prerank), value = prerank),
        path = path, col.names = FALSE
    )
}

# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
write_gmt <- function(gs, path) {
    data <- as.data.frame(gs)
    data$genesets <- vapply(data$genesets, paste, character(1L),
        collapse = "\t", USE.NAMES = FALSE
    )
    write_table(data, path, col.names = FALSE)
}

read_gmt <- function(path, ...) {
    data <- strsplit(read_lines(path), "\t", fixed = TRUE)
    parse_genesets(path, unclass(data))
}

# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29
read_gmx <- function(path, ...) {
    data <- read_table(path, header = FALSE, sep = "\t")
    parse_genesets(path, unclass(data))
}

parse_genesets <- function(path, genesets) {
    genesets <- genesets[lengths(genesets) > 2L]
    if (length(genesets) == 0L) cli::cli_abort("No genesets in {.path {path}}")
    ids <- vapply(genesets, .subset, 1L, character(1L), USE.NAMES = FALSE)
    descriptions <- vapply(
        genesets, .subset, 2L, character(1L),
        USE.NAMES = FALSE
    )
    genesets <- lapply(genesets, function(geneset) {
        trimws(.subset(geneset, -(1:2)), "both")
    })
    new_genesets(genesets, ids = ids, descriptions = descriptions)
}
