##' Build virus metadata for seqcombo
##'
##' Create a `virus_info` data frame for plotting reassortment diagrams.
##'
##' @param id unique virus identifiers
##' @param segment_color a list-column of segment colors, one entry per virus
##' @param segment_name optional list-column of segment names, one entry per virus
##' @param x,y optional x/y coordinates
##' @param virus_size optional relative virus sizes
##' @param label optional text labels
##' @param label_position optional label positions; one of `left`, `right`,
##'   `below`, `upper`, or `none`
##' @param ... additional columns to attach to the returned data frame
##' @return a `virus_info` data frame
##' @export
##' @author Guangchuang Yu
build_virus_info <- function(id, segment_color, segment_name = NULL, x = NULL, y = NULL,
                             virus_size = NULL, label = NULL,
                             label_position = NULL, ...) {
    n <- length(id)
    segment_color <- normalize_segment_color(segment_color, n)
    segment_name <- normalize_segment_name(segment_name, n)

    virus_info <- data.frame(
        id = id,
        x = normalize_optional_column(x, n, NA_real_, "x"),
        y = normalize_optional_column(y, n, NA_real_, "y"),
        segment_color = I(segment_color),
        virus_size = normalize_optional_column(virus_size, n, 1, "virus_size"),
        stringsAsFactors = FALSE
    )
    if (!is.null(segment_name)) {
        virus_info$segment_name <- I(segment_name)
    }

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
    validate_segment_name(virus_info)
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


##' Build reassortment edge metadata from long-format records
##'
##' Aggregate long-format reassortment records into a `flow_info` data frame.
##'
##' @param data a data frame with one row per reassortment record
##' @param from column name containing source virus identifiers
##' @param to column name containing target virus identifiers
##' @param segment optional column name containing segment identifiers
##' @param weight optional column name containing edge weights to sum
##' @param keep optional character vector of additional per-edge columns to keep
##' @return a `flow_info` data frame
##' @export
##' @author Guangchuang Yu
build_flow_info_from_long <- function(data, from = "from", to = "to",
                                      segment = NULL, weight = NULL,
                                      keep = NULL) {
    require_columns(data, c(from, to))

    key <- paste(data[[from]], data[[to]], sep = "\r")
    split_data <- split(data, key, drop = TRUE)

    edge_info <- data.frame(
        from = vapply(split_data, function(d) d[[from]][1], FUN.VALUE = data[[from]][1]),
        to = vapply(split_data, function(d) d[[to]][1], FUN.VALUE = data[[to]][1]),
        stringsAsFactors = FALSE
    )

    args <- list(from = edge_info$from, to = edge_info$to)

    if (!is.null(segment)) {
        require_columns(data, segment)
        args$segment_name <- lapply(split_data, function(d) unique(as.character(d[[segment]])))
    }

    if (!is.null(weight)) {
        require_columns(data, weight)
        args$weight <- vapply(split_data, function(d) sum(d[[weight]], na.rm = TRUE), FUN.VALUE = numeric(1))
    }

    if (!is.null(keep) && length(keep) > 0) {
        require_columns(data, keep)
        for (col in keep) {
            args[[col]] <- vapply(split_data, function(d) unique_scalar(d[[col]], col), FUN.VALUE = data[[col]][1])
        }
    }

    do.call(build_flow_info, args)
}


##' Build virus metadata from long-format segment records
##'
##' Aggregate long-format segment records into a `virus_info` data frame.
##'
##' @param data a data frame with one row per virus-segment pair
##' @param id column name containing virus identifiers
##' @param segment column name containing segment identifiers
##' @param color column name containing segment colors
##' @param x,y optional column names for coordinates
##' @param virus_size optional column name for relative virus sizes
##' @param label optional column name for text labels
##' @param label_position optional column name for label positions
##' @param keep optional character vector of additional per-virus columns to keep
##' @param segment_order optional vector specifying segment order
##' @return a `virus_info` data frame
##' @export
##' @author Guangchuang Yu
build_virus_info_from_long <- function(data, id = "id", segment = "segment",
                                       color = "color", x = NULL, y = NULL,
                                       virus_size = NULL, label = NULL,
                                       label_position = NULL, keep = NULL,
                                       segment_order = NULL) {
    require_columns(data, c(id, segment, color))

    if (is.null(segment_order)) {
        segment_order <- unique(data[[segment]])
    }

    order_index <- match(data[[segment]], segment_order)
    if (any(is.na(order_index))) {
        stop("all segment values must be present in 'segment_order'...")
    }

    data <- data[order(data[[id]], order_index), , drop = FALSE]
    split_data <- split(data, data[[id]], drop = TRUE)

    segment_color <- lapply(split_data, function(d) d[[color]])
    segment_name <- lapply(split_data, function(d) as.character(d[[segment]]))
    ids <- names(split_data)

    per_virus <- data.frame(id = ids, stringsAsFactors = FALSE)
    scalar_fields <- c(
        x = x,
        y = y,
        virus_size = virus_size,
        label = label,
        label_position = label_position
    )

    for (nm in names(scalar_fields)) {
        col <- scalar_fields[[nm]]
        if (!is.null(col)) {
            require_columns(data, col)
            per_virus[[nm]] <- vapply(split_data, function(d) unique_scalar(d[[col]], col), FUN.VALUE = data[[col]][1])
        }
    }

    if (!is.null(keep) && length(keep) > 0) {
        require_columns(data, keep)
        for (col in keep) {
            per_virus[[col]] <- vapply(split_data, function(d) unique_scalar(d[[col]], col), FUN.VALUE = data[[col]][1])
        }
    }

    args <- list(
        id = per_virus$id,
        segment_color = segment_color,
        segment_name = segment_name,
        x = per_virus$x,
        y = per_virus$y,
        virus_size = per_virus$virus_size,
        label = per_virus$label,
        label_position = per_virus$label_position
    )

    extra_names <- setdiff(
        colnames(per_virus),
        c("id", "x", "y", "virus_size", "label", "label_position")
    )
    if (length(extra_names) > 0) {
        args <- c(args, as.list(per_virus[extra_names]))
    }

    do.call(build_virus_info, args)
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


normalize_segment_name <- function(segment_name, n) {
    if (is.null(segment_name)) {
        return(NULL)
    }
    if (!is.list(segment_name)) {
        if (length(segment_name) != n) {
            stop("'segment_name' must be a list with one entry per virus...")
        }
        segment_name <- as.list(segment_name)
    }

    if (length(segment_name) != n) {
        stop("'segment_name' must be a list with one entry per virus...")
    }

    lapply(segment_name, as.character)
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


require_columns <- function(data, columns) {
    missing_col <- setdiff(columns, colnames(data))
    if (length(missing_col) > 0) {
        stop(sprintf("missing required columns: %s", paste(missing_col, collapse = ", ")))
    }
}


unique_scalar <- function(x, name) {
    x <- unique(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(NA)
    }
    if (length(x) > 1) {
        stop(sprintf("column '%s' must have a single value per virus...", name))
    }
    x[[1]]
}


##' Build a segment caption string
##'
##' Create a compact caption that lists segment names in order.
##'
##' @param segment_name a character vector of segment names
##' @param prefix prefix text used before the segment names
##' @return a single caption string
##' @export
##' @author Guangchuang Yu
build_segment_caption <- function(segment_name, prefix = "Gene segments: ") {
    if (is.list(segment_name)) {
        segment_name <- segment_name[[1]]
    }
    segment_name <- as.character(segment_name)
    if (length(segment_name) == 0) {
        stop("'segment_name' must contain at least one segment name...")
    }
    paste0(prefix, paste(segment_name, collapse = ", "))
}
