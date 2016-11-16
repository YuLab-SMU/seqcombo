##' @rdname plot-methods
##' @exportMethod plot
##' @param width bin width
##' @param title plot title
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom grid unit.pmax
##' @importFrom cowplot plot_grid
##' @author guangchuang yu
setMethod("plot", signature(x="SeqDiff"),
          function(x, width=50, title="auto") {
              nn <- names(x@sequence)
              if (is.null(title) || is.na(title)) {
                  title <- ""
              } else if (title == "auto") {
                  title <- paste(nn[-x@reference], "nucelotide differences relative to", nn[x@reference])
              }

              p1 <- plot_difference_count(x@diff, width) + ggtitle(title)
              p2 <- plot_difference(x@diff)

              gp1<- ggplot_gtable(ggplot_build(p1))
              gp2<- ggplot_gtable(ggplot_build(p2))
              maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
              gp1$widths[2:3] <- maxWidth
              gp2$widths[2:3] <- maxWidth

              ## pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.7, 0.7), "null"))))
              ## gp2$vp = viewport(layout.pos.row = 2, layout.pos.col = 1)
              ## grid.draw(gp2)
              ## gp1$vp = viewport(layout.pos.row = 1, layout.pos.col = 1)
              ## grid.draw(gp1)
              plot_grid(gp1, gp2, ncol=1)
          }
          )



##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 theme_minimal
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
plot_difference <- function(x) {
    yy = 4:1
    names(yy) = c("A", "C", "G", "T")
    x$y <- yy[x$difference]
    ggplot(x, aes_(x=~position, y=~y, color=~difference)) +
        geom_segment(aes_(x=~position, xend=~position, y=~y, yend=~y+.8)) +
        xlab("Nucleotide Position") + ylab(NULL) +
        scale_y_continuous(breaks=yy, labels=names(yy)) +
        theme_minimal() +
        theme(legend.position="none")+
        theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
}

##' @importFrom ggplot2 geom_col
##' @importFrom ggplot2 theme_bw
plot_difference_count <- function(x, width) {
    d <- nucleotide_difference_count(x, width)
    ggplot(d, aes_(x=~position, y=~count)) + geom_col(fill='red', width=width*.8) +
        xlab(NULL) + ylab("Difference") + theme_bw()
}

