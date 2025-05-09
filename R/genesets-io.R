# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
write_gmt <- function(gs, path) {
    data <- as.data.frame(gs)[c("ids", "descriptions", "genesets")]
    data$descriptions <- replace_na(data$descriptions, "")
    lines <- vapply(vec_seq_along(data), function(i) {
        paste(unlist(vec_slice(data, i), TRUE, FALSE), collapse = "\t")
    }, character(1L), USE.NAMES = FALSE)
    write_lines(lines, path)
}

read_gmt <- function(path, ...) {
    data <- strsplit(read_lines(path), "\t", fixed = TRUE)
    parse_genesets(path, data)
}

# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29
read_gmx <- function(path, ...) {
    data <- read_table(path, header = FALSE, sep = "\t")
    parse_genesets(path, as.list(data))
}

parse_genesets <- function(path, genesets) {
    genesets <- genesets[lengths(genesets) > 2L]
    if (length(genesets) == 0L) cli::cli_abort("No genesets in {.path {path}}")
    ids <- vapply(genesets, .subset, character(1L), 1L, USE.NAMES = FALSE)
    descriptions <- vapply(
        genesets, .subset, character(1L), 2L,
        USE.NAMES = FALSE
    )
    descriptions[!nzchar(descriptions) | descriptions == "NA"] <- NA_character_
    genesets <- lapply(genesets, function(geneset) {
        trimws(.subset(geneset, -(1:2)), "both")
    })
    new_genesets(genesets, ids = ids, descriptions = descriptions)
}
