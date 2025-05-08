#' GSEA prerank algorithm
#'
#' @name GSEAPrerank
NULL

#' @describeIn GSEAPrerank R-based preranked GSEA algorithm (official). See
#' [`GSEA`][GSEA::GSEA] for details.
GSEAPrerank <- new_class("GSEAPrerank",
    properties = list(
        exponent = prop_number_decimal(
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.double(value)
                prop(self, "exponent") <- value
                self
            },
            default = 1,
            allow_infinite = FALSE
        ),
        nperm = prop_number_whole(
            min = 1,
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.integer(value)
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

# GSEA method (official) -------------------------------
method(bridge, GSEAPrerank) <- function(source, target, method) {
    check_remote_installed(
        "GSEA", "GSEA-MSigDB/GSEA_R",
        "to use {.field GSEAPrerank} method"
    )

    # Prepare the output directory
    odir <- tempdir()
    dir_create(odir)

    # Prepare the input file
    rnk <- tempfile("rnk", fileext = ".rnk")
    write_rnk(source, rnk)
    on.exit(file.remove(rnk), add = TRUE)
    gmt <- tempfile("gmt", fileext = ".gmt")
    write_gmt(target, gmt)
    on.exit(file.remove(gmt), add = TRUE)

    # run GSEA
    getExportedValue("GSEA", "GSEA")(
        input.ds = rnk,
        gs.db = gmt,
        output.directory = odir,
        reshuffling.type = "gene.labels",
        nperm = method@nperm,
        weighted.score.type = method@exponent,
        fdr.q.val.threshold = 1,
        adjust.FDR.q.val = TRUE,
        gs.size.threshold.min = 0L,
        gs.size.threshold.max = length(source),
        preproc.type = 0L,
        random.seed = sample.int(1e6L, 1L),
        save.intermediate.results = FALSE,
        use.fast.enrichment.routine = TRUE
    )
    odir
}

# prerank method ---------------------------------------
#' @describeIn GSEAPrerank Rust-based preranked GSEA algorithm.
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
#' @include utils-S7.R
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
