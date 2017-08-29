## ' access diff slot
## '
## ' @param x SeqDiff object
## ' @return data.frame
## ' @examples
## ' fas <- list.files(system.file("examples","GVariation", package="seqcombo"),
## ' pattern="fas", full.names=TRUE)
## ' x1 <- seqdiff(fas[1], reference=1)
## ' diff(x1)
##' @method diff SeqDiff
##' @export
diff.SeqDiff <- function(x, ...) {
    x@diff
}


