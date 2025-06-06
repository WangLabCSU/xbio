% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genesets-gs.R
\name{gs_ids}
\alias{gs_ids}
\alias{gs_terms}
\alias{gs_descs}
\title{Extract attributes from a geneset or genesets object}
\usage{
gs_ids(gs)

gs_terms(gs)

gs_descs(gs)
}
\arguments{
\item{gs}{A \code{\link[=geneset]{geneset()}} or \code{\link[=genesets]{genesets()}} object.}
}
\value{
A character vector containing the extracted attributes.
If \code{gs} is a \code{geneset} object, a character scalar is returned.
}
\description{
These functions extract specific attributes from \code{\link[=geneset]{geneset()}} or
\code{\link[=genesets]{genesets()}} objects:
\itemize{
\item IDs (\code{gs_ids()})
\item terms (\code{gs_terms()})
\item descriptions (\code{gs_descs()}).
}
}
