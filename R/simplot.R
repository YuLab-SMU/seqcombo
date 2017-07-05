##' Sequence similarity plot
##'
##'
##' @title simplot
##' @param file alignment fast file
##' @param query query sequence
##' @return ggplot object
##' @importFrom Biostrings readDNAStringSet
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 ggtitle
##' @export
##' @author guangchuang yu
##' @examples
##' fas <- system.file("examples/GVariation/sample_alignment.fa", package="seqcombo")
##' simplot(fas, 'CF_YL21')
simplot <- function(file, query) {
    aln <- readDNAStringSet(file)
    nn <- names(aln)
    idx <- which(nn != query)
    res <- do.call(rbind, lapply(idx, function(i) {
        x <- toCharacter(aln[i]) == toCharacter(aln[query])
        pos <- seq_along(x)
        y <- data.frame(sequence=nn[i], position=pos, similarity=cumsum(x)/pos * 100)
    }))
    ggplot(res, aes_(x=~position, y=~similarity)) +
        geom_line(aes_(group=~sequence, color=~sequence)) +
        xlab("Nucleotide Position") + ylab("Similarity (%)") +
        ggtitle(paste("Sequence similarities compare to", query)) +
        theme_minimal() +
        theme(legend.title=element_blank(), legend.position=c(.95, .15))
}


toCharacter <- function(x) {
	unlist(strsplit(toString(x),""))
}


