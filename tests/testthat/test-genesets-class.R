test_that("creates new_genesets with named gene sets", {
    gs <- list(a = c("X", "Y"), b = c("Y", "Z"))
    obj <- new_genesets(gs)
    expect_s3_class(obj, "enricher_genesets")
    expect_named(obj, c("a", "b"))
})

test_that("errors when terms are missing and genesets are unnamed", {
    gs <- list(c("A", "B"), c("C", "D"))
    expect_error(
        new_genesets(gs),
        regexp = "must be provided or"
    )
})

test_that("assigns terms as names when provided", {
    gs <- list(c("A", "B"), c("C", "D"))
    terms <- c("Set1", "Set2")
    obj <- new_genesets(gs, terms = terms)
    expect_named(obj, terms)
})

test_that("errors when terms and genesets have mismatched lengths", {
    gs <- list(c("A", "B"), c("C", "D"))
    terms <- c("Set1")
    expect_error(
        new_genesets(gs, terms = terms),
        regexp = "must have the same length"
    )
})

test_that("assigns descriptions when provided", {
    gs <- list(c("A", "B"), c("C", "D"))
    terms <- c("Set1", "Set2")
    desc <- c("desc1", "desc2")
    obj <- new_genesets(gs, terms = terms, descriptions = desc)
    expect_s3_class(obj, "enricher_genesets")
    expect_named(obj, terms)
    expect_identical(attr(obj, "descriptions"), desc)
})

test_that("errors when descriptions and genesets have mismatched lengths", {
    gs <- list(c("A", "B"), c("C", "D"))
    terms <- c("Set1", "Set2")
    desc <- c("desc1") # too short
    expect_error(
        new_genesets(gs, terms = terms, descriptions = desc),
        regexp = "must have the same length"
    )
})

# Name Handling
test_that("names<-.enricher_genesets cworks for valid input", {
    gs <- list(c("A", "B"), c("C"))
    obj <- new_genesets(gs,
        terms = c("old1", "old2"),
        descriptions = c("d1", "d2")
    )

    names(obj) <- c("new1", "new2")
    expect_named(obj, c("new1", "new2"))
})

test_that("names<-.enricher_genesets preserves descriptions when renaming", {
    gs <- list(c("A", "B"), c("C"))
    obj <- new_genesets(gs,
        terms = c("a", "b"),
        descriptions = c("desc1", "desc2")
    )

    names(obj) <- c("x", "y")
    expect_named(obj, c("x", "y"))
    expect_equal(attr(obj, "descriptions"), c("desc1", "desc2"))
})

test_that("names<-.enricher_genesets rejects NULL in names", {
    gs <- list(a = c("A", "B"))
    obj <- new_genesets(gs)
    expect_error(
        names(obj) <- NULL,
        "Cannot remove the names of genesets"
    )
})

test_that("names<-.enricher_genesets rejects NA in names", {
    x <- new_genesets(list(a = c("TP53", "EGFR"), b = c("BRCA1", "BRCA2")))
    expect_error(names(x) <- c("Valid", NA), "Names cannot be missing or empty")
})

test_that("names<-.enricher_genesets rejects empty string in names", {
    x <- new_genesets(list(a = c("TP53", "EGFR"), b = c("BRCA1", "BRCA2")))
    new_names <- c("", "Valid")
    expect_error(names(x) <- new_names, "Names cannot be missing or empty")
})

# Subsetting and indexing
test_that("`[` subset returns enriched object", {
    gs <- list(a = c("A", "B"), b = c("C"))
    desc <- c("Pathway A", "Pathway B")
    obj <- new_genesets(gs, descriptions = desc)
    sub <- obj[1]

    expect_s3_class(sub, "enricher_genesets")
    expect_equal(length(sub), 1)
    expect_equal(names(sub), "a")
    expect_equal(attr(sub, "descriptions"), "Pathway A")
})

test_that("`[[` and `$` retain description", {
    gs <- list(a = c("A", "B"))
    obj <- new_genesets(gs, descriptions = "desc")
    expect_equal(attr(obj[["a"]], "description"), "desc")
    expect_equal(attr(obj$a, "description"), "desc")
})

# Type conversion and casting
test_that("conversion to data.frame works", {
    gs <- list(a = c("A", "B"), b = c("C"))
    desc <- c("Pathway A", "Pathway B")
    obj <- new_genesets(gs, descriptions = desc)

    df <- as.data.frame(obj)
    expect_s3_class(df, "data.frame")
    expect_equal(ncol(df), 3)
    expect_named(df, c("terms", "descriptions", "genesets"))
})

test_that("conversion to list preserves names", {
    gs <- list(a = c("X"), b = c("Y", "Z"))
    desc <- c("one", "two")
    obj <- new_genesets(gs, descriptions = desc)
    lst <- as.list(obj)
    expect_type(lst, "list")
    expect_named(lst, names(obj))
})

# rep() and length<-
test_that("rep replicates gene sets correctly", {
    gs <- list(a = c("A", "B"), b = "C")
    obj <- new_genesets(gs, descriptions = c("d1", "d2"))

    r <- rep(obj, 2)
    expect_equal(length(r), 4L)
    expect_s3_class(r, "enricher_genesets")
})

test_that("length<- truncates the object", {
    gs <- list(a = c("A", "B"), b = "C", c = "D")
    obj <- new_genesets(gs, descriptions = c("d1", "d2", "d3"))
    length(obj) <- 2
    expect_equal(names(obj), c("a", "b"))
    expect_equal(attr(obj, "descriptions"), c("d1", "d2"))
})

test_that("length<- errors if increasing length", {
    gs <- list(a = c("A", "B"))
    obj <- new_genesets(gs)
    expect_error(length(obj) <- 2, "Cannot set length greater")
})
