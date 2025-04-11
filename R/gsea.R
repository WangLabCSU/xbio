gsea <- function(object, gs, ..., method, min_size = NULL, max_size = NULL) {
    dots <- rlang::dots_list(...,
        .ignore_empty = "all", .named = NULL,
        .homonyms = "error"
    )
    if (missing(method) || is.null(method)) {
        if (S7_inherits(.subset2(dots, 1L), gseaMethod)) {
            method <- .subset2(dots, 1L)
            dots <- dots[-1L]
        } else {
            method <- gseaMultilevel()
        }
    }
    if (!rlang::is_named2(dots)) {
        cli::cli_abort("All input in {.arg ...} must be named")
    }
    assert_number_whole(min_size, min = 1, allow_null = TRUE)
    assert_number_whole(min_size, min = 1, allow_null = TRUE)
    object <- prerank(object)
    gs <- genesets(gs)
    gs <- filter_genesets(gs, min_size, max_size)
    run_gsea(method, object = object, gs = gs, params = dots)
}

prerank <- function(x, ...) UseMethod("prerank")

run_gsea <- new_generic(
    "run_gsea", "method",
    function(method, object, gs, params) S7_dispatch()
)

gseaMethod <- new_class("gseaMethod")

#' @include utils-S7.R
gseaSimple <- new_class(
    "gseaSimple", gseaMethod,
    properties = list(
        nperm = prop_number_whole(
            min = 1,
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.integer(value)
                prop(self, "nperm") <- value
                self
            },
            default = 1L
        ),
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
        ),
        exponential = prop_number_decimal(
            setter = function(self, value) {
                if (is.numeric(value)) value <- as.double(value)
                prop(self, "exponential") <- value
                self
            },
            default = 1,
            allow_infinite = FALSE
        )
    )
)

method(run_gsea, gseaSimple) <- function(method, object, gs, params) {
    if (length(params) > 0L) props(method) <- params
    fgsea::fgseaSimple(
        gs,
        object,
        minSize = 1,
        maxSize = length(object),
        scoreType = method@score_type,
        nperm = method@nperm,
        nproc = method@threads,
        gseaParam = method@exponential
    )
}

gseaMultilevel <- new_class(
    "gseaMultilevel",
    gseaSimple,
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
    )
)

method(run_gsea, gseaMultilevel) <- function(method, object, gs, params) {
    if (length(params) > 0L) props(method) <- params
    fgsea::fgseaMultilevel(
        gs,
        object,
        minSize = 1,
        maxSize = length(object),
        sampleSize = method@sample_size,
        eps = method@eps,
        scoreType = method@score_type,
        nPermSimple = method@nperm,
        nproc = method@threads,
        gseaParam = method@exponential
    )
}
