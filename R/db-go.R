godb <- function(database, strategy = NULL, save = NULL, verbose = TRUE) {

}

godb_download <- function(database) {
    if (identical(database, "ontology")) {
        url <- "https://current.geneontology.org/ontology/go-basic.obo"
    } else {
        url <- sprintf(
            "http://current.geneontology.org/annotations/%s.gaf.gz",
            database
        )
    }
    ofile <- file.path(cache_dir(), "GO", basename(url))
    utils::download.file(url, ofile)
    if (identical(database, "ontology")) {
        read_lines(ofile)
    } else {
        read_table(ofile)
    }
}

GO_ONTOLOGY <- c("BP", "CC", "MF")
