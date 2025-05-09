#' GSEA prerank algorithm
#'
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

S4_register(GSEAPrerank)

method(repr_source, GSEAPrerank) <- function(method, source) {
    repr_metrics(source, `_arg` = "source")
}

method(repr_target, GSEAPrerank) <- function(method, target) {
    repr_genesets(target, `_arg` = "target")
}

# prerank method ---------------------------------------
#' @describeIn GSEAPrerank Rust-based GSEA gene-permutation algorithm.
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

S4_register(GSEAGene)

method(bridge, GSEAGene) <- function(source, target, method) {
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
    indices <- out$indices
    out$indices <- NULL
    new_data_frame(c(
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

S4_register(GSEASimple)

method(bridge, GSEASimple) <- function(source, target, method) {
    check_bioc_installed("fgsea", "to use {.field GSEASimple} method")
    fgsea::fgseaSimple(
        target,
        source,
        minSize = 0,
        maxSize = length(source),
        scoreType = method@score_type,
        nperm = method@nperm,
        nproc = method@threads,
        gseaParam = method@exponent
    )
}

#' @describeIn GSEAPrerank `fgsea` algorithm. See
#' [`fgseaMultilevel`][fgsea::fgseaMultilevel] for details.
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

S4_register(GSEAMultilevel)

method(bridge, GSEAMultilevel) <- function(source, target, method) {
    check_bioc_installed("fgsea", "to use {.field GSEAMultilevel} method")
    fgsea::fgseaMultilevel(
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
}

# GSEA method (official) -------------------------------
#' @describeIn GSEAPrerank R-based GSEA gene-permutation algorithm (official).
#' See [`GSEA`][GSEA::GSEA] for details.
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

S4_register(GSEABroadGene)

method(bridge, GSEABroadGene) <- function(source, target, method) {
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
    out <- rename(out, c(
        ES = "es", NES = "nes",
        `NOM p-val` = "pvalue", `FDR q-val` = "fdr"
    ))
    out <- vec_cbind(vec_slice(gs_data, index), out)
    vec_slice(out, order(index))
}
