##' Set layout for reassortment plot
##'
##' @title set_layout
##' @param virus_info virus information
##' @param flow_info flow information
##' @param layout layout method
##' @param preserve_x whether to preserve existing x coordinates and only update y
##' @param preserve_y whether to preserve existing y coordinates and only update x
##' @param spread whether to evenly spread viruses sharing the same preserved
##' coordinate along the other axis (useful when many viruses share the same
##' time point)
##' @return updated `virus_info`
##' @importFrom igraph graph_from_data_frame
##' @importFrom igraph V
##' @importFrom igraph layout_nicely
##' @importFrom yulab.utils get_fun_from_pkg
##' @export
##' @author Guangchuang Yu
set_layout <- function(virus_info, flow_info, layout = "layout_nicely",
                       preserve_x = FALSE, preserve_y = FALSE,
                       spread = FALSE) {
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

    g <- graph_from_data_frame(flow_info[, c("from", "to")])
    coord <- layout(g)
    i <- match(as.character(V(g)$name), virus_info$id)

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
    if (spread && preserve_x) {
        virus_info$y[i] <- spread_grouped_values(layout_y, virus_info$x[i])
    }
    if (spread && preserve_y) {
        virus_info$x[i] <- spread_grouped_values(layout_x, virus_info$y[i])
    }

    virus_info
}


##' Arrange values into evenly spaced slots, grouped by `groups`
##'
##' Nodes sharing a group occupy consecutive slots (with a gap between
##' groups), preserving the relative order of the original values. Useful to
##' spread viruses that share the same time point on a timeline layout.
##'
##' @noRd
spread_grouped_values <- function(values, groups, gap = 1) {
    n <- length(values)
    stopifnot(n == length(groups))

    groups_chr <- as.character(groups)
    groups_chr[is.na(groups_chr)] <- "__NA__"
    levels_sorted <- sort_unique_key(groups_chr)
    ord <- order(match(groups_chr, levels_sorted), values)

    breaks <- rle(groups_chr[ord])
    res <- numeric(n)
    cursor <- 0
    start <- 1
    for (k in breaks$lengths) {
        idx <- ord[start:(start + k - 1)]
        res[idx] <- cursor + seq_len(k)
        cursor <- cursor + k + gap
        start <- start + k
    }
    ## flip so that the first group appears at the top after plotting
    max_res <- max(res)
    res <- max_res + 1 - res
    res
}

sort_unique_key <- function(chr) {
    u <- unique(chr)
    num <- suppressWarnings(as.numeric(u))
    if (!anyNA(num)) {
        return(u[order(num)])
    }
    sort(u)
}


##' Layout reassortment plots on a timeline
##'
##' Keep temporal coordinates on one axis and automatically place viruses on the
##' other axis.
##'
##' @title layout_timeline
##' @param virus_info virus information
##' @param flow_info flow information
##' @param time_col column name in `virus_info` that contains temporal order
##' @param axis which axis should preserve time, one of `x` or `y`
##' @param decreasing whether to reverse the automatically computed axis
##' @param spread whether to evenly spread viruses sharing the same time point
##' along the other axis
##' @return updated `virus_info`
##' @export
##' @author Guangchuang Yu
layout_timeline <- function(virus_info, flow_info, time_col = "x",
                            axis = c("x", "y"), decreasing = FALSE,
                            spread = TRUE) {
    axis <- match.arg(axis)
    if (!time_col %in% colnames(virus_info)) {
        stop(sprintf("'%s' column is required in 'virus_info'...", time_col))
    }

    preserve_x <- axis == "x"
    preserve_y <- axis == "y"
    virus_info <- set_layout(
        virus_info = virus_info,
        flow_info = flow_info,
        preserve_x = preserve_x,
        preserve_y = preserve_y,
        spread = spread
    )

    virus_info[[axis]] <- virus_info[[time_col]]
    other_axis <- if (axis == "x") "y" else "x"
    if (decreasing) {
        virus_info[[other_axis]] <- max(virus_info[[other_axis]], na.rm = TRUE) -
            virus_info[[other_axis]]
    }

    virus_info
}
