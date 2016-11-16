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
setMethod("show",signature(object="SeqDiff"),
          function(object) {
              cat("sequence differences of", paste0(names(object@sequence), collapse=" and "), '\n')
              d <- object@diff$difference %>% table %>% as.data.frame
              cat(sum(d$Freq), "sites differ:\n")
              freq <- d[,2]
              names(freq) <- d[,1]
              print(freq)
          })
