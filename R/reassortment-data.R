##' Build virus metadata for seqcombo
##'
##' Create a `virus_info` data frame for plotting reassortment diagrams.
##'
##' @param id unique virus identifiers
##' @param segment_color a list-column of segment colors, one entry per virus
##' @param x,y optional x/y coordinates
##' @param virus_size optional relative virus sizes
##' @param label optional text labels
##' @param label_position optional label positions; one of `left`, `right`,
##'   `below`, `upper`, or `none`
##' @param ... additional columns to attach to the returned data frame
##' @return a `virus_info` data frame
##' @export
##' @author Guangchuang Yu
build_virus_info <- function(id, segment_color, x = NULL, y = NULL,
                             virus_size = NULL, label = NULL,
                             label_position = NULL, ...) {
    n <- length(id)
    segment_color <- normalize_segment_color(segment_color, n)

    virus_info <- data.frame(
        id = id,
        x = normalize_optional_column(x, n, NA_real_, "x"),
        y = normalize_optional_column(y, n, NA_real_, "y"),
        segment_color = I(segment_color),
        virus_size = normalize_optional_column(virus_size, n, 1, "virus_size"),
        stringsAsFactors = FALSE
    )

    if (!is.null(label) || !is.null(label_position)) {
        virus_info$label <- normalize_optional_column(label, n, NA_character_, "label")
        virus_info$label_position <- normalize_label_position(label_position, n)
    }

    extra_columns <- list(...)
    if (length(extra_columns) > 0) {
        for (nm in names(extra_columns)) {
            virus_info[[nm]] <- normalize_optional_column(extra_columns[[nm]], n, NULL, nm)
        }
    }

    validate_virus_info(virus_info, require_coordinates = FALSE)
    validate_segment_color(virus_info)
    virus_info
}


##' Build reassortment edge metadata for seqcombo
##'
##' Create a `flow_info` data frame for plotting reassortment flows.
##'
##' @param from source virus identifiers
##' @param to target virus identifiers
##' @param ... additional columns to attach to the returned data frame
##' @return a `flow_info` data frame
##' @export
##' @author Guangchuang Yu
build_flow_info <- function(from, to, ...) {
    n <- length(from)
    if (length(to) != n) {
        stop("'from' and 'to' must have the same length...")
    }

    flow_info <- data.frame(from = from, to = to, stringsAsFactors = FALSE)

    extra_columns <- list(...)
    if (length(extra_columns) > 0) {
        for (nm in names(extra_columns)) {
            flow_info[[nm]] <- normalize_optional_column(extra_columns[[nm]], n, NULL, nm)
        }
    }

    validate_flow_info(flow_info, virus_id = NULL)
    flow_info
}


normalize_optional_column <- function(value, n, default = NULL, name) {
    if (is.null(value)) {
        if (is.null(default)) {
            return(NULL)
        }
        if (length(default) == 1L) {
            return(rep(default, n))
        }
        return(default)
    }

    if (length(value) == 1L) {
        return(rep(value, n))
    }

    if (length(value) != n) {
        stop(sprintf("'%s' must have length 1 or %d...", name, n))
    }

    value
}


normalize_segment_color <- function(segment_color, n) {
    if (!is.list(segment_color)) {
        if (length(segment_color) != n) {
            stop("'segment_color' must be a list with one entry per virus...")
        }
        segment_color <- as.list(segment_color)
    }

    if (length(segment_color) != n) {
        stop("'segment_color' must be a list with one entry per virus...")
    }

    lapply(segment_color, as.character)
}


normalize_label_position <- function(label_position, n) {
    label_position <- normalize_optional_column(label_position, n, "none", "label_position")
    label_position <- tolower(label_position)

    supported <- c("left", "right", "below", "upper", "none")
    invalid <- setdiff(unique(label_position), supported)
    if (length(invalid) > 0) {
        stop("'label_position' must be one of 'left', 'right', 'below', 'upper', or 'none'...")
    }

    label_position
}
