#' Gene Set Enrichment Analysis (GSEA) - Preranked Variants
#'
#' @description
#' A family of enrichment methods based on the GSEA preranked algorithm,
#' evaluating gene set enrichment using a ranked gene list as input. These
#' methods differ in implementation (e.g., R, Rust, or C++ backends), p-value
#' estimation techniques (e.g., empirical, adaptive multilevel), and performance
#' optimizations.
#'
#' @param exponent A positive numeric scalar that controls the weighting of the
#' ranking metric. A value of `1` corresponds to the original GSEA method.
#'
#' @param nperm A positive integer specifying the number of permutations to
#' estimate statistical significance. Higher values improve precision but
#' increase runtime.
#' @name GSEAPrerank
NULL

#' @include utils-S7.R
GSEAPrerank <- new_class("GSEAPrerank",
    properties = list(
        exponent = prop_number_decimal(
            setter = function(self, value) {
                if (is.integer(value)) value <- as.double(value)
                prop(self, "exponent") <- value
                self
            },
            default = 1,
            allow_infinite = FALSE
        ),
        nperm = prop_number_whole(
            min = 1,
            setter = function(self, value) {
                if (is.double(value)) value <- as.integer(value)
                prop(self, "nperm") <- value
                self
            },
            default = 1000L
        )
    )
)

method(repr_source, list(GSEAPrerank, class_any)) <- function(method, source) {
    repr_metrics(source, `_arg` = "source")
}

method(repr_target, list(GSEAPrerank, class_any)) <- function(method, target) {
    repr_genesets(target, `_arg` = "target")
}

# prerank method ---------------------------------------
# Adaptive Tail Approximation for GSEA (ATA-GSEA)
# The standard empirical P-value in GSEA is slow and inaccurate for small
# P-values because:
# 1. It relies on many permutations (often 1,000–100,000),
# 2. For small observed enrichment scores, extreme tail behavior is poorly
#    sampled
# 3. Memory and runtime become limiting for large datasets.
#' @describeIn GSEAPrerank Rust-based GSEA gene-permutation algorithm.
#' @param threads Integer. Number of threads to use for parallel computation.
#'   Defaults to all available cores.
#' @export
GSEAGene <- new_class("GSEAGene", GSEAPrerank,
    properties = list(
        threads = prop_number_whole(
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.integer(value)
                prop(self, "threads") <- value
                self
            },
            min = 1,
            max = quote(as.double(parallel::detectCores()))
        )
    ),
    constructor = function(exponent = 1, nperm = 1000, threads = NULL) {
        new_object(
            GSEAPrerank(exponent = exponent, nperm = nperm),
            threads = threads %||% parallel::detectCores()
        )
    }
)

method(bridge, GSEAGene) <- function(method, source, target) {
    out <- rust_call(
        "gsea_gene_permutate",
        names(source),
        source,
        genesets = c(unclass(target)), # remove attributes
        exponent = method@exponent,
        nperm = method@nperm,
        threads = method@threads,
        seed = sample.int(1e6L, 1L)
    )
    new_gsea_result(c(
        list(
            ids = names(target),
            terms = gs_terms(target),
            descriptions = gs_descs(target)
        ),
        out
    ))
}

# fgsea method -----------------------------------------
#' @describeIn GSEAPrerank `fgsea` algorithm. See
#' [`fgseaSimple`][fgsea::fgseaSimple] for details.
#' @param score_type Character string indicating how to handle gene-level
#'   statistic signs when calculating enrichment scores. One of:
#'   \itemize{
#'     \item `"std"`: Use original signed statistics (standard behavior).
#'     \item `"pos"`: Use only positive statistics; negative values set to zero.
#'     \item `"neg"`: Use only negative statistics; positive values set to zero.
#'   }
#' @export
GSEASimple <- new_class(
    "GSEASimple", GSEAPrerank,
    properties = list(
        score_type = prop_match(c("std", "pos", "neg")),
        threads = prop_number_whole(
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.integer(value)
                prop(self, "threads") <- value
                self
            },
            min = 1,
            max = quote(as.double(parallel::detectCores())),
            default = 1L
        )
    )
)

method(bridge, GSEASimple) <- function(method, source, target) {
    check_bioc_installed("fgsea", "to use {.field GSEASimple} method")
    # Save the genesets information, will be restored in the result
    gs_data <- new_data_frame(list(
        ids = names(target),
        terms = gs_terms(target),
        descriptions = gs_descs(target)
    ))
    names(target) <- as.character(seq_len(vec_size(target)))

    # Run GSEA
    out <- fgsea::fgseaSimple(
        target,
        source,
        minSize = 0,
        maxSize = length(source),
        scoreType = method@score_type,
        nperm = method@nperm,
        nproc = method@threads,
        gseaParam = method@exponent
    )

    # restore the original ordering
    out <- as.data.frame(out)
    index <- as.integer(out$pathway)
    out <- out[c("ES", "NES", "pval", "padj", "leadingEdge")]
    names(out) <- c("es", "nes", "pvalue", "fdr", "leading_edge")
    out <- vec_cbind(vec_slice(gs_data, index), out)
    new_gsea_result(vec_slice(out, order(index)))
}

#' @describeIn GSEAPrerank `fgsea` algorithm. See
#' [`fgseaMultilevel`][fgsea::fgseaMultilevel] for details.
#' @param sample_size Number of samples per level for multilevel splitting.
#' @param eps Minimum p-value to estimate; must be a small positive number.
#' @export
GSEAMultilevel <- new_class(
    "GSEAMultilevel",
    GSEASimple,
    properties = list(
        sample_size = prop_number_whole(
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.integer(value)
                prop(self, "sample_size") <- value
                self
            },
            min = 1,
            default = 101L
        ),
        eps = prop_number_decimal(
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.double(value)
                prop(self, "eps") <- value
                self
            },
            default = 1e-50,
            allow_infinite = FALSE
        )
    ),
    constructor = function(eps = 1e-50, exponent = 1, nperm = 1000L,
                           sample_size = 101L, score_type = "std",
                           threads = 1L) {
        new_object(
            GSEASimple(
                nperm = nperm, score_type = score_type,
                threads = threads, exponent = exponent
            ),
            sample_size = sample_size,
            eps = eps
        )
    }
)

method(bridge, GSEAMultilevel) <- function(method, source, target) {
    check_bioc_installed("fgsea", "to use {.field GSEAMultilevel} method")
    # Save the genesets information, will be restored in the result
    gs_data <- new_data_frame(list(
        ids = names(target),
        terms = gs_terms(target),
        descriptions = gs_descs(target)
    ))
    names(target) <- as.character(seq_len(vec_size(target)))

    # Run GSEA
    out <- fgsea::fgseaMultilevel(
        target,
        source,
        minSize = 0,
        maxSize = length(source),
        sampleSize = method@sample_size,
        eps = method@eps,
        scoreType = method@score_type,
        nPermSimple = method@nperm,
        nproc = method@threads,
        gseaParam = method@exponent
    )

    # restore the original ordering
    out <- as.data.frame(out)
    index <- as.integer(out$pathway)
    out <- out[c("ES", "NES", "pval", "padj", "leadingEdge")]
    names(out) <- c("es", "nes", "pvalue", "fdr", "leading_edge")
    out <- vec_cbind(vec_slice(gs_data, index), out)
    new_gsea_result(vec_slice(out, order(index)))
}

# GSEA method (official) -------------------------------
#' @describeIn GSEAPrerank R-based GSEA gene-permutation algorithm (official).
#' See [`GSEA`][GSEA::GSEA] for details.
#' @param odir Character string giving the output directory where the
#'   GSEA results will be saved. If `NA` (default), a temporary directory
#'   is used and deleted automatically after processing. If a path is
#'   provided, the directory will be created if it does not exist and
#'   retained after completion.
#' @export
GSEABroadGene <- new_class(
    "GSEABroadGene", GSEAPrerank,
    properties = list(
        odir = prop_string(
            allow_empty = FALSE, allow_na = TRUE,
            default = NA_character_
        )
    )
)

method(bridge, GSEABroadGene) <- function(method, source, target) {
    check_remote_installed(
        "GSEA", "GSEA-MSigDB/GSEA_R",
        "to use {.field GSEABroadGene} method"
    )

    # Prepare the output directory
    if (is.null(odir <- method@odir) || is.na(odir)) {
        odir <- tempfile("GSEABroadGene")
        dir_create(odir)
        on.exit(unlink(odir, recursive = TRUE, force = TRUE), add = TRUE)
    } else {
        dir_create(odir)
    }

    # Save the genesets information, will be restored in the result
    gs_data <- new_data_frame(list(
        ids = names(target),
        terms = gs_terms(target),
        descriptions = gs_descs(target)
    ))

    # Prepare the input file
    rnk <- tempfile("rnk", fileext = ".rnk")
    write_rnk(source, rnk)
    on.exit(file.remove(rnk), add = TRUE)
    names(target) <- as.character(seq_len(vec_size(target)))
    gmt <- tempfile("gmt", fileext = ".gmt")
    write_gmt(target, gmt)
    on.exit(file.remove(gmt), add = TRUE)

    # run GSEA
    out <- getExportedValue("GSEA", "GSEA")(
        input.ds = rnk,
        input.cls = list(),
        gs.db = gmt,
        output.directory = odir,
        gsea.type = "preranked",
        reshuffling.type = "gene.labels",
        nperm = method@nperm,
        weighted.score.type = method@exponent,
        fdr.q.val.threshold = 1,
        adjust.FDR.q.val = FALSE,
        gs.size.threshold.min = 0L,
        gs.size.threshold.max = length(source),
        preproc.type = 0L,
        random.seed = sample.int(1e6L, 1L),
        save.intermediate.results = FALSE,
        use.fast.enrichment.routine = TRUE
    )
    # result <- list.files(odir, pattern = "NA_pos\\.txt", full.names = TRUE)
    # out <- read_table(result, comment = "", header = TRUE)
    out <- vec_rbind(!!!out)
    index <- as.integer(out$GS)
    out <- out[c("ES", "NES", "NOM p-val", "FDR q-val")]
    names(out) <- c("es", "nes", "pvalue", "fdr")
    out <- vec_cbind(vec_slice(gs_data, index), out)
    new_gsea_result(vec_slice(out, order(index)))
}
