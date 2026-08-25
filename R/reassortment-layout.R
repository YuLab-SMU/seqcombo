##' Set layout for reassortment plot
##'
##' @title set_layout
##' @param virus_info virus information
##' @param flow_info flow information
##' @param layout layout method
##' @param preserve_x whether to preserve existing x coordinates and only update y
##' @param preserve_y whether to preserve existing y coordinates and only update x
##' @return updated `virus_info`
##' @importFrom igraph graph.data.frame
##' @importFrom igraph V
##' @importFrom igraph layout.auto
##' @importFrom yulab.utils get_fun_from_pkg
##' @export
##' @author Guangchuang Yu
set_layout <- function(virus_info, flow_info, layout = "layout.auto",
                       preserve_x = FALSE, preserve_y = FALSE) {
    if (preserve_x && preserve_y) {
        stop("'preserve_x' and 'preserve_y' cannot both be TRUE...")
    }

    validate_virus_info(virus_info, require_coordinates = FALSE)
    validate_flow_info(flow_info, virus_info$id)

    if (preserve_x && !"x" %in% colnames(virus_info)) {
        stop("'x' column is required when 'preserve_x = TRUE'...")
    }
    if (preserve_y && !"y" %in% colnames(virus_info)) {
        stop("'y' column is required when 'preserve_y = TRUE'...")
    }

    if (is.character(layout)) {
        layout <- get_fun_from_pkg("igraph", layout)
    }

    g <- graph.data.frame(flow_info[, c("from", "to")])
    coord <- layout(g)
    i <- match(as.character(V(g)), virus_info$id)

    if (!"x" %in% colnames(virus_info)) {
        virus_info$x <- NA_real_
    }
    if (!"y" %in% colnames(virus_info)) {
        virus_info$y <- NA_real_
    }

    layout_x <- max(coord[, 1]) - coord[, 1]
    layout_y <- max(coord[, 2]) - coord[, 2]

    if (!preserve_x) {
        virus_info$x[i] <- layout_x
    }
    if (!preserve_y) {
        virus_info$y[i] <- layout_y
    }

    virus_info
}
