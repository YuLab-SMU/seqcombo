
##' plot method generics
##'
##'
##' @docType methods
##' @name plot
##' @rdname plot-methods
##' @title plot method
##' @param x object
##' @param ... Additional argument list
##' @return plot
##' @importFrom stats4 plot
##' @export
if ( !isGeneric("plot") )
    setGeneric("plot", function(x, ...) standardGeneric("plot"))
