gs_map <- function(gs, key_source, key_target,
                   annodb = NULL, organism = NULL) {
    assert_s3_class(gs, "enricher_genesets")
    assert_string(key_source, allow_empty = FALSE)
    assert_string(key_target, allow_empty = FALSE)

    if (vec_unique_count(attr(gs, "organism")) > 1L) {
        cli::cli_abort(c(
            "Mapping genesets from multiple {.arg organism} values is not supported.",
            i = "Please map genesets for each organism separately before merging."
        ))
    }

    if (vec_size(gs) == 0L) return(gs) # styler: off

    # Infer the annodb ----------------------------
    if (is.null(annodb)) {
        if (is.null(organism)) {
            organism <- vec_unique(attr(gs, "organism"))
            if (is.null(organism)) {
                cli::cli_abort(paste(
                    "Please specify {.arg annodb} to map {.arg gs} into",
                    "{.field {key_target}}, or specify {.arg organism}",
                    "so we can infer {.arg annodb} for you."
                ))
            }
        } else if (!.rlang_check_string(organism, allow_empty = FALSE)) {
            cli::cli_abort("{.arg organism} must be a single string")
        }

        annodb <- switch(organism,
            hsa = "org.Hs.eg.db",
            cli::cli_abort(c(
                "Cannot infer {.arg annodb} for {.field {organism}}",
                i = "Please provide {.arg annodb} manually"
            ))
        )
    }

    if (is.character(annodb)) {
        check_bioc_installed(annodb)
        annodb <- getExportedValue(annodb, annodb)
    }

    # mapping the the genes in genesets into keytype
    gs_lapply(gs, function(geneset) {
        if (length(geneset) == 0L) return(geneset) # styler: off
        out <- AnnotationDbi::mapIds(
            annodb,
            keys = geneset,
            column = key_target,
            keytype = key_source
        )
        out[is.na(out) | out == ""] <- NA_character_
        out
    })
}

#' @keywords internal
#' @noRd
gs_trim <- function(gs, ...) {
    gs_lapply(gs, function(geneset) geneset[!is.na(geneset) & geneset != ""])
}

gs_filter <- function(gs, min_size = NULL, max_size = NULL) {
    assert_s3_class(gs, "enricher_genesets")
    assert_number_whole(min_size, min = 1, allow_null = TRUE)
    assert_number_whole(min_size, min = 1, allow_null = TRUE)
    if (is.null(min_size) && is.null(max_size)) {
        return(gs)
    }
    sizes <- list_sizes(gs)
    if (!is.null(min_size) && !is.null(max_size)) {
        keep <- sizes >= min_size & sizes <= max_size
        out_pattern <- sprintf("[%d, %d]", min_size, max_size)
    } else if (!is.null(min_size)) {
        keep <- sizes >= min_size
        out_pattern <- sprintf("[%d, Inf)", min_size)
    } else {
        keep <- sizes <= max_size
        out_pattern <- sprintf("(0, %d]", max_size)
    }
    if (!all(keep)) {
        cli::cli_warn(sprintf(
            "Removing {sum(!keep)} gene set{?s} out of %s",
            out_pattern
        ))
        gs <- gs[keep]
    }
    gs
}

gs_lapply <- function(gs, ...) {
    data <- vec_proxy(gs)
    data$genesets <- lapply(data$genesets, ...)
    vec_restore(data, gs)
}
