% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-show.R
\docType{methods}
\name{show}
\alias{show}
\alias{SeqDiff-class}
\alias{show,SeqDiff-method}
\title{show method}
\usage{
show(object)
}
\arguments{
\item{object}{SeqDiff object}
}
\value{
message
}
\description{
show method
}
\examples{
fas <- list.files(system.file("examples","GVariation", package="seqcombo"), pattern="fas", full.names=TRUE)
x1 <- seqdiff(fas[1], reference=1)
x1
}
