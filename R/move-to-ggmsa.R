setClass("SeqDiff",
         representation = representation(
             file = "character",
             sequence = "BStringSet",
             reference = "numeric",
             diff = "data.frame"
         )
         )


##' @method diff SeqDiff
##' @export
diff.SeqDiff <- function(x, ...) {
    x@diff
}


##' @export
##' @importFrom ggplot2 qplot
##' @examples
##' fas <- list.files(system.file("examples","GVariation", package="seqcombo"),
##'                  pattern="fas", full.names=TRUE)
##' x <- lapply(fas, seqdiff)
##' plts <- lapply(x, plot)
##' plot_grid(plotlist=plts, ncol=1, labels=LETTERS[1:3])
cowplot::plot_grid


##' plot method for SeqDiff object
##'
##' @name plot
##' @rdname plot-methods
##' @exportMethod plot
##' @aliases plot,SeqDiff,ANY-method
##' @docType methods
##' @param x SeqDiff object
##' @param width bin width
##' @param title plot title
##' @param xlab xlab
##' @param by one of 'bar' and 'area'
##' @param fill fill color of upper part of the plot
##' @param colors color of lower part of the plot
##' @param xlim limits of x-axis
##' @return plot
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom grid unit.pmax
##' @importFrom cowplot plot_grid
##' @author guangchuang yu
##' @examples
##' fas <- list.files(system.file("examples","GVariation", package="seqcombo"), pattern="fas", full.names=TRUE)
##' x1 <- seqdiff(fas[1], reference=1)
##' plot(x1)
setMethod("plot", signature(x="SeqDiff"),
          function(x, width=50, title="auto",
                   xlab = "Nucleotide Position",
                   by="bar", fill="firebrick",
                   colors=c(A="#E495A5", C="#ABB065", G="#39BEB1", T="#ACA4E2"),
                   xlim = NULL) {
              message("This function is deprecated and will be removed in next release...")
              nn <- names(x@sequence)
              if (is.null(title) || is.na(title)) {
                  title <- ""
              } else if (title == "auto") {
                  title <- paste(nn[-x@reference], "nucelotide differences relative to", nn[x@reference])
              }

              p1 <- plot_difference_count(x@diff, width, by=by, fill=fill) + ggtitle(title)
              p2 <- plot_difference(x@diff, colors=colors, xlab)

              if (!is.null(xlim)) {
                  p1 <- p1 + xlim(xlim)
                  p2 <- p2 + xlim(xlim)
              }

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
              plot_grid(gp1, gp2, ncol=1, rel_heights=c(.7, .4))
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
##' @importFrom ggplot2 scale_color_manual
plot_difference <- function(x, colors, xlab="Nucleotide Position") {
    yy = 4:1
    names(yy) = c("A", "C", "G", "T")
    x$y <- yy[x$difference]
    n <- sum(is.na(x$y))
    if (n > 0) {
        message(n, " sites contain deletions or ambiguous bases, which will be ignored in current implementation...")
    }
    x <- x[!is.na(x$y),]
    p <- ggplot(x, aes_(x=~position, y=~y, color=~difference))

    p + geom_segment(aes_(x=~position, xend=~position, y=~y, yend=~y+.8)) +
        xlab(xlab) + ylab(NULL) +
        scale_y_continuous(breaks=yy, labels=names(yy)) +
        theme_minimal() +
        theme(legend.position="none")+
        theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
        scale_color_manual(values=colors)
}

##' @importFrom ggplot2 geom_col
##' @importFrom ggplot2 geom_area
##' @importFrom ggplot2 theme_bw
plot_difference_count <- function(x, width, by = 'bar', fill='red') {
    by <- match.arg(by, c("bar", "area"))
    if (by == 'bar') {
        geom <- geom_col(fill=fill, width=width)
        keep0 <- FALSE
    } else if (by == "area") {
        geom <- geom_area(fill=fill)
        keep0 <- TRUE
    }
    d <- nucleotide_difference_count(x, width, keep0)
    p <- ggplot(d, aes_(x=~position, y=~count))
    p + geom + xlab(NULL) + ylab("Difference") + theme_bw()
}


##' show method
##'
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##' @title show method
##' @param object SeqDiff object
##' @return message
##' @importFrom methods show
##' @exportMethod show
##' @aliases SeqDiff-class
##'   show,SeqDiff-method
##' @usage show(object)
##' @examples
##' fas <- list.files(system.file("examples","GVariation", package="seqcombo"), pattern="fas", full.names=TRUE)
##' x1 <- seqdiff(fas[1], reference=1)
##' x1
setMethod("show",signature(object="SeqDiff"),
          function(object) {
              message("This function is deprecated and will be removed in next release...")
              cat("sequence differences of", paste0(names(object@sequence), collapse=" and "), '\n')
              d <- object@diff$difference %>% table %>% as.data.frame
              cat(sum(d$Freq), "sites differ:\n")
              freq <- d[,2]
              names(freq) <- d[,1]
              print(freq)
          })



##' calculate difference of two aligned sequences
##'
##'
##' @title seqdiff
##' @param fasta fasta file
##' @param reference which sequence serve as reference, 1 or 2
##' @return SeqDiff object
##' @export
##' @importFrom Biostrings readBStringSet
##' @importClassesFrom Biostrings BStringSet
##' @importFrom methods new
##' @author guangchuang yu
##' @examples
##' fas <- list.files(system.file("examples","GVariation", package="seqcombo"), pattern="fas", full.names=TRUE)
##' seqdiff(fas[1], reference=1)
seqdiff <- function(fasta, reference=1) {
    message("This function is deprecated and will be removed in next release...")

    sequence <- readBStringSet(fasta)
    if (length(sequence) != 2 && length(width(sequence)) != 1) {
        stop("fas should contains 2 aligned sequences...")
    }
    diff <- nucleotide_difference(sequence, reference)
    new("SeqDiff",
        file = fasta,
        sequence = sequence,
        reference = reference,
        diff = diff)
}

##' @importFrom magrittr %>%
##' @importFrom Biostrings toString
##' @importFrom Biostrings width
nucleotide_difference <- function(x, reference=1) {
    n <- width(x[1])
    nn <- seq_len(n)
    s1 <- x[1] %>% toString %>% substring(nn, nn)
    s2 <- x[2] %>% toString %>% substring(nn, nn)

    pos <- which(s1 != s2)
    if (reference == 1) {
        diff <- s2[pos]
    } else {
        diff <- s1[pos]
    }

    return(data.frame(position = pos,
                      difference = diff,
                      stringsAsFactors = FALSE))
}




##' @importFrom dplyr group_by
##' @importFrom dplyr summarize
##' @importFrom dplyr select
##' @importFrom dplyr n
nucleotide_difference_count <- function(x, width=50, keep0=FALSE) {
    n <- max(x$position)
    bin <- rep(1:ceiling(n/width), each=width)
    position <- c(seq_len(n)[!duplicated(bin)], n)
    x$bin <- bin[x$pos]
    y <- x %>% group_by(bin) %>%
        summarize(position=min(position), count = n()) %>%
        select(-bin)
    y$position <- position[findInterval(y$position, position)]
    if (keep0) {
        itv <- seq(1, n, width)
        yy <- data.frame(position = itv[!itv %in% y$position],
                         count = 0)
        y <- rbind(y, yy)
        y <- y[order(y$position, decreasing=FALSE),]
    }
    return(y)
}



##' Sequence similarity plot
##'
##'
##' @title simplot
##' @param file alignment fast file
##' @param query query sequence
##' @param window sliding window size (bp)
##' @param step step size to slide the window (bp)
##' @param group whether grouping sequence
##' @param id position to extract id for grouping; only works if group = TRUE
##' @param sep separator to split sequence name; only works if group = TRUE
##' @param sd whether display standard deviation of similarity among each group; only works if group=TRUE
##' @return ggplot object
##' @importFrom Biostrings readDNAStringSet
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 geom_ribbon
##' @importFrom magrittr %<>%
##' @importFrom dplyr group_by_
##' @importFrom dplyr summarize_
##' @export
##' @author guangchuang yu
##' @examples
##' fas <- system.file("examples/GVariation/sample_alignment.fa", package="seqcombo")
##' simplot(fas, 'CF_YL21')
simplot <- function(file, query, window=200, step=20, group=FALSE, id, sep, sd=FALSE) {
    message("This function is deprecated and will be removed in next release...")
    aln <- readDNAStringSet(file)
    nn <- names(aln)
    if (group) {
        g <- vapply(strsplit(nn, sep), function(x) x[id], character(1))
    }

    idx <- which(nn != query)
    w <- width(aln[query])
    start <- seq(1, w, by=step)
    end <- start + window - 1
    start <- start[end <= w]
    end <- end[end <= w]
    res <- lapply(idx, function(i) {
        x <- toCharacter(aln[i]) == toCharacter(aln[query])
        ## pos <- seq_along(x)
        ## data.frame(sequence=nn[i], position=pos, similarity=cumsum(x)/pos * 100)
        pos <- round((start+end)/2)
        sim <- vapply(seq_along(start), function(j) {
            mean(x[start[j]:end[j]])
        }, numeric(1))
        ## sim <- c(cummean(x[1:(pos[1]-1)]), sim, cummean(x[(pos[length(pos)]+1):w]))
        ## pos <- c(1:(pos[1]-1), pos, (pos[length(pos)]+1):w)

        y <- data.frame(sequence=nn[i], position = pos, similarity = sim)
        if(group) {
            y$group <- g[i]
        }
        return(y)
    }) %>% do.call(rbind, .)

    if (group) {
        res %<>% group_by_(~position, ~group) %>%
            summarize_(msim=~mean(similarity), sd=~sd(similarity))
    }


    if (group) {
        p <- ggplot(res, aes_(x=~position, y=~msim, group=~group))
        if (sd) p <- p + geom_ribbon(aes_(ymin=~msim-sd, ymax=~msim+sd, fill=~group), alpha=.25) #fill='grey70')
        p <- p + geom_line(aes_(color=~group))
    } else {
        p <- ggplot(res, aes_(x=~position, y=~similarity)) +
            geom_line(aes_(group=~sequence, color=~sequence))
    }

    p + xlab("Nucleotide Position") + ylab("Similarity (%)") +
        ggtitle(paste("Sequence similarities compare to", query)) +
        theme_minimal() +
        theme(legend.title=element_blank()) #, legend.position=c(.95, .15))
}


toCharacter <- function(x) {
	unlist(strsplit(toString(x),""))
}


