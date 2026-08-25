validate_virus_info <- function(virus_info, require_coordinates = TRUE) {
    require_col <- c("id", "segment_color")
    missing_col <- setdiff(require_col, colnames(virus_info))
    if (length(missing_col) > 0) {
        stop("'id' and 'segment_color' columns are required in 'virus_info'...")
    }

    if (require_coordinates) {
        coord_col <- c("x", "y")
        missing_coord <- setdiff(coord_col, colnames(virus_info))
        if (length(missing_coord) > 0) {
            stop("'x' and 'y' columns are required in 'virus_info'...")
        }
    }
}


validate_flow_info <- function(flow_info, virus_id = NULL) {
    if (!all(c("from", "to") %in% colnames(flow_info))) {
        stop("'from' and 'to' columns are required in 'flow_info'...")
    }

    if (is.null(virus_id)) {
        return(invisible(TRUE))
    }

    flow_id <- unique(c(flow_info$from, flow_info$to))
    missing_id <- setdiff(flow_id, virus_id)
    if (length(missing_id) > 0) {
        stop("all 'from' and 'to' ids in 'flow_info' must exist in 'virus_info$id'...")
    }
}


validate_segment_color <- function(virus_info) {
    segment_length <- lengths(virus_info$segment_color)
    if (any(segment_length == 0)) {
        stop("each 'segment_color' entry in 'virus_info' must contain at least one color...")
    }
}


validate_link_style <- function(link_style) {
    match.arg(link_style, c("segment", "curve"))
}
