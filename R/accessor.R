##' access diff slot
##'
##' @param x SeqDiff object
##' @return data.frame
##' @examples
##' fas <- list.files(system.file("examples","GVariation", package="seqcombo"),
##' pattern="fas", full.names=TRUE)
##' x1 <- seqdiff(fas[1], reference=1)
##' get_diff(x1)
##' @export
get_diff <- function(x) {
    UseMethod("get_diff", x)
}


##' @method get_diff SeqDiff
##' @export
get_diff.SeqDiff <- function(x) x@diff

