gsea <- new_generic("gsea", c("object", "genesets"))

# Add a method for our class
method(gsea, list(class_any, class_list)) <- function(object, genesets) {
    object
}
