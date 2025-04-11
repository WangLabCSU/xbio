gsea <- new_generic(
    "gsea", "method",
    function(method, object, gs, ..., min_size = NULL, max_size = NULL) {
        object <- prerank(object)
        gs <- genesets(gs)
        gs <- filter_genesets(gs, min_size, max_size)
        S7_dispatch()
    }
)

prerank <- function(x, ...) UseMethod("prerank")

#' @include utils-S7.R
gseaSimple <- new_class(
    "gseaSimple",
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
        scoreType = prop_match(c("std", "pos", "neg")),
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

method(gsea, gseaSimple) <- function(method, object, gs, ...,
                                     min_size = NULL, max_size = NULL) {
    fgsea::fgseaSimple(
        gs,
        object,
        minSize = 1,
        maxSize = length(object),
        scoreType = method@method,
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

method(gsea, gseaMultilevel) <- function(method, object, gs, ...,
                                         min_size = NULL, max_size = NULL) {
    fgsea::fgseaMultilevel(
        gs,
        object,
        minSize = 1,
        maxSize = length(object),
        sampleSize = method@sample_size,
        eps = method@eps,
        scoreType = method@method,
        nPermSimple = method@nperm,
        nproc = method@threads,
        gseaParam = method@exponential
    )
}
