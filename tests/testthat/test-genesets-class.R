test_that("creates new_genesets with named gene sets", {
    gs <- list(a = c("X", "Y"), b = c("Y", "Z"))
    obj <- new_genesets(gs)
    expect_s3_class(obj, "xbio_genesets")
    expect_named(obj, c("a", "b"))
})

test_that("assigns ids as names when provided", {
    gs <- list(c("A", "B"), c("C", "D"))
    ids <- c("Set1", "Set2")
    obj <- new_genesets(gs, ids = ids)
    expect_named(obj, ids)
})

test_that("errors when ids and genesets have mismatched lengths", {
    gs <- list(c("A", "B"), c("C", "D"))
    ids <- c("Set1")
    expect_error(
        new_genesets(gs, ids = ids),
        regexp = "must have the same length"
    )
})

test_that("assigns descriptions when provided", {
    gs <- list(c("A", "B"), c("C", "D"))
    ids <- c("Set1", "Set2")
    desc <- c("desc1", "desc2")
    obj <- new_genesets(gs, ids = ids, descriptions = desc)
    expect_s3_class(obj, "xbio_genesets")
    expect_named(obj, ids)
    expect_identical(gs_descs(obj), desc)
})

test_that("errors when descriptions and genesets have mismatched lengths", {
    gs <- list(c("A", "B"), c("C", "D"))
    ids <- c("Set1", "Set2")
    desc <- c("desc1") # too short
    expect_error(
        new_genesets(gs, ids = ids, descriptions = desc),
        regexp = "must have the same length"
    )
})

# Name Handling
test_that("names<-.xbio_genesets cworks for valid input", {
    gs <- list(c("A", "B"), c("C"))
    obj <- new_genesets(gs,
        ids = c("old1", "old2"),
        descriptions = c("d1", "d2")
    )

    names(obj) <- c("new1", "new2")
    expect_named(obj, c("new1", "new2"))
})

test_that("names<-.xbio_genesets preserves descriptions when renaming", {
    gs <- list(c("A", "B"), c("C"))
    obj <- new_genesets(gs,
        ids = c("a", "b"),
        descriptions = c("desc1", "desc2")
    )

    names(obj) <- c("x", "y")
    expect_named(obj, c("x", "y"))
    expect_equal(gs_descs(obj), c("desc1", "desc2"))
})

# Subsetting and indexing
test_that("`[` subset returns enriched object", {
    gs <- list(a = c("A", "B"), b = c("C"))
    desc <- c("Pathway A", "Pathway B")
    obj <- new_genesets(gs, descriptions = desc)
    sub <- obj[1]

    expect_s3_class(sub, "xbio_genesets")
    expect_equal(length(sub), 1)
    expect_equal(names(sub), "a")
    expect_equal(gs_descs(sub), "Pathway A")
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
    expect_equal(ncol(df), 4L)
    expect_named(df, c("ids", "terms", "descriptions", "genesets"))
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
    expect_s3_class(r, "xbio_genesets")
})

test_that("length<- truncates the object", {
    gs <- list(a = c("A", "B"), b = "C", c = "D")
    obj <- new_genesets(gs, descriptions = c("d1", "d2", "d3"))
    length(obj) <- 2
    expect_equal(names(obj), c("a", "b"))
    expect_equal(gs_descs(obj), c("d1", "d2"))
})

test_that("length<- errors if increasing length", {
    gs <- list(a = c("A", "B"))
    obj <- new_genesets(gs)
    expect_error(length(obj) <- 2, "Cannot set length greater")
})
