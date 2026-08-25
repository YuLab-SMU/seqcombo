##' @importFrom ggplot2 geom_polygon
##' @importFrom ggplot2 geom_rect
##' @importFrom ggplot2 aes
##' @importFrom rlang .data
##' @importFrom utils modifyList
geom_virus_capsule <- function(mapping, virus_info, hex_data, color, fill,
                               alpha = 0.5, linewidth = 1) {
    hex.df <- do.call("rbind", hex_data)

    if (typeof(color) == "language") {
        vcol <- all.vars(color)
        if (!vcol %in% colnames(virus_info)) {
            stop("color variable not available...")
        }
        hex.df[, vcol] <- virus_info[[vcol]][match(hex.df$id, virus_info$id)]
        mapping <- modifyList(mapping, aes(color = .data[[vcol]]))
    }

    if (typeof(fill) == "language") {
        vf <- all.vars(fill)
        if (!vf %in% colnames(virus_info)) {
            stop("fill variable not available...")
        }
        hex.df[, vf] <- virus_info[[vf]][match(hex.df$id, virus_info$id)]
        mapping <- modifyList(mapping, aes(fill = .data[[vf]]))
    }
    mapping <- modifyList(mapping, aes(group = .data[["id"]]))

    params <- list(
        mapping = mapping,
        data = hex.df,
        alpha = alpha,
        linewidth = linewidth,
        inherit.aes = FALSE
    )
    if (typeof(color) == "character") {
        params <- modifyList(list(color = color), params)
    }
    if (typeof(fill) == "character") {
        params <- modifyList(list(fill = fill), params)
    }

    do.call(geom_polygon, params)
}


geom_gene_segment <- function(hexd, color, height = 0.68, g_height = 0.65,
                              g_width = 0.8) {
    n <- length(color)
    y <- hexd$y
    yh <- diff(range(y)) / 4 * g_height / 0.5
    y <- y[y >= mean(y) - yh & y <= mean(y) + yh]

    ymin <- ymax <- seq(min(y), max(y), length.out = n + 1)
    ymin <- ymin[-(n + 1)]
    ymax <- ymax[-1]
    adjust <- (ymax - ymin) * (1 - height) / 2
    ymin <- ymin + adjust
    ymax <- ymax - adjust

    x <- hexd$x
    xx <- x[hexd$y == max(y)]
    xmin <- min(x)
    xmax <- max(x)
    xadj <- (xmax - xmin) * (1 - g_width) / 2

    d <- data.frame(
        xmin = max(xmin + xadj, min(xx)),
        xmax = min(xmax - xadj, max(xx)),
        ymin = ymin,
        ymax = ymax,
        color = rev(color)
    )

    dd <- split(d, d$color)
    lapply(seq_along(dd), function(i) {
        geom_rect(
            aes(xmin = .data[["xmin"]], ymin = .data[["ymin"]],
                xmax = .data[["xmax"]], ymax = .data[["ymax"]]),
            data = dd[[i]],
            fill = dd[[i]]$color[1],
            inherit.aes = FALSE,
            show.legend = FALSE
        )
    })
}


set_virus_size <- function(virus_info, ASP, v_shape = "ellipse") {
    v_shape <- match.arg(v_shape, c("hexagon", "ellipse"))

    if (!"virus_size" %in% colnames(virus_info)) {
        virus_info$virus_size <- 1
    }

    if (ASP < 1) {
        virus_info$virus_size <- virus_info$virus_size / 20 * diff(range(virus_info$y, na.rm = TRUE))
    } else {
        virus_info$virus_size <- virus_info$virus_size / 20 * diff(range(virus_info$x, na.rm = TRUE))
    }

    if (v_shape == "ellipse") {
        virus_info$virus_size <- virus_info$virus_size * 0.5
    }

    virus_info
}


get_capsule_data <- function(virus_info, ASP, v_shape) {
    lapply(seq_len(nrow(virus_info)), function(i) {
        d <- generate_capsule_data(
            x = virus_info$x[i],
            y = virus_info$y[i],
            size = virus_info$virus_size[i],
            ASP = ASP,
            shape = v_shape
        )
        d$id <- virus_info$id[i]
        d
    })
}


generate_capsule_data <- function(x, y, size, ASP, shape) {
    if (shape == "ellipse") {
        return(generate_ellipse_data(x, y, size, ASP))
    }
    generate_hex_data(x, y, size, ASP)
}


generate_hex_data <- function(x, y, size, ASP) {
    asp <- estimate_asp(ASP)
    data.frame(
        x = c(rep(-sqrt(3) / 2, 2), 0, rep(sqrt(3) / 2, 2), 0) * size * asp[1] + x,
        y = c(0.5, -0.5, -1, -0.5, 0.5, 1) * size * asp[2] + y
    )
}


generate_ellipse_data <- function(x, y, size, ASP) {
    asp <- estimate_asp(ASP)
    a <- 3
    b <- 4
    xx <- seq(-sqrt(a), sqrt(a), length.out = 500)
    yy <- sqrt(b * (1 - xx^2 / a))
    data.frame(
        x = c(xx, rev(xx)) * size * asp[1] + x,
        y = c(yy, rev(-yy)) * size * asp[2] + y
    )
}


route_reassortment_edges <- function(virus_info, flow_info, hex_data, ASP = 1) {
    x <- virus_info$x
    y <- virus_info$y

    x.from <- x[match(flow_info$from, virus_info$id)]
    x.to <- x[match(flow_info$to, virus_info$id)]
    y.from <- y[match(flow_info$from, virus_info$id)]
    y.to <- y[match(flow_info$to, virus_info$id)]

    width <- vapply(hex_data, function(d) max(d$x), numeric(1)) - x
    height <- vapply(hex_data, function(d) max(d$y), numeric(1)) - y
    names(width) <- names(height) <- virus_info$id

    xadj <- diff(range(virus_info$x)) * 0.01
    yadj <- diff(range(virus_info$y)) * 0.01

    xdiff <- (x.to - x.from) / xadj
    ydiff <- (y.to - y.from) / yadj
    vertical_idx <- abs(ydiff) > abs(xdiff)

    x.from.adj <- width[flow_info$from] + xadj
    x.from.adj[vertical_idx] <- 0
    x.to.adj <- width[flow_info$to] + xadj
    x.to.adj[vertical_idx] <- 0

    y.from.adj <- rep(0, nrow(flow_info))
    y.from.adj[vertical_idx] <- height[flow_info$from][vertical_idx] + yadj
    y.to.adj <- rep(0, nrow(flow_info))
    y.to.adj[vertical_idx] <- height[flow_info$to][vertical_idx] + yadj

    x.direction <- sign(x.to - x.from)
    x.from <- x.from + x.from.adj * x.direction
    x.to <- x.to - x.to.adj * x.direction

    y.direction <- sign(y.to - y.from)
    y.from <- y.from + y.from.adj * y.direction
    y.to <- y.to - y.to.adj * y.direction

    edge_data <- data.frame(x = x.from, xend = x.to, y = y.from, yend = y.to)

    edge_data <- adjust_edge_endpoints(edge_data, "x", "y", width, height, flow_info$from, vertical_idx)
    edge_data <- adjust_edge_endpoints(edge_data, "xend", "yend", width, height, flow_info$to, vertical_idx)
    edge_data <- adjust_shared_targets(edge_data, flow_info, width, height, vertical_idx)

    edge_data[, c("x", "y", "xend", "yend")]
}


adjust_edge_endpoints <- function(edge_data, x_col, y_col, width, height, id, vertical_idx) {
    edge_data$id <- id
    dup <- which(duplicated(edge_data[, c(x_col, y_col, "id")]))
    if (length(dup) == 0) {
        edge_data$id <- NULL
        return(edge_data)
    }

    dup <- dup[!duplicated(edge_data[dup, "id"])]
    for (i in dup) {
        j <- edge_data$id == edge_data$id[i] &
            edge_data[, x_col] == edge_data[i, x_col] &
            edge_data[, y_col] == edge_data[i, y_col]

        if (all(vertical_idx[j])) {
            w <- width[edge_data$id[j]][1] / 2
            offset <- seq(-w, w, length.out = sum(j) + 2)
            offset <- offset[-c(1, length(offset))]
            if (x_col == "x") {
                order_idx <- order(edge_data[j, "xend"], decreasing = FALSE)
            } else {
                order_idx <- order(edge_data[j, "x"], decreasing = FALSE)
            }
            edge_data[j, x_col] <- edge_data[j, x_col] + offset[order_idx]
        } else {
            h <- height[edge_data$id[j]][1] / 2
            offset <- seq(-h, h, length.out = sum(j) + 2)
            offset <- offset[-c(1, length(offset))]
            if (y_col == "y") {
                order_idx <- order(edge_data[j, "yend"], decreasing = FALSE)
            } else {
                order_idx <- order(edge_data[j, "y"], decreasing = FALSE)
            }
            edge_data[j, y_col] <- edge_data[j, y_col] + offset[order_idx]
        }
    }

    edge_data$id <- NULL
    edge_data
}


adjust_shared_targets <- function(edge_data, flow_info, width, height, vertical_idx) {
    edge_data$from <- flow_info$from
    edge_data$to <- flow_info$to
    paired_idx <- match(edge_data$to, edge_data$from)
    shared_idx <- which(edge_data$xend == edge_data$x[paired_idx] &
                        edge_data$yend == edge_data$y[paired_idx])

    if (length(shared_idx) > 0) {
        shared_idx <- shared_idx[!duplicated(edge_data$to[shared_idx])]

        for (j in shared_idx) {
            if (vertical_idx[j]) {
                h <- height[edge_data$to[j]] / 2
                offset <- seq(-h, h, length.out = 4)
                if (edge_data$x[j] > edge_data$xend[paired_idx[j]]) {
                    edge_data$xend[j] <- edge_data$xend[j] + offset[3]
                    edge_data$x[paired_idx[j]] <- edge_data$x[paired_idx[j]] + offset[2]
                } else {
                    edge_data$xend[j] <- edge_data$xend[j] + offset[2]
                    edge_data$x[paired_idx[j]] <- edge_data$x[paired_idx[j]] + offset[3]
                }
            } else {
                w <- width[edge_data$to[j]] / 2
                offset <- seq(-w, w, length.out = 4)
                if (edge_data$y[j] > edge_data$yend[paired_idx[j]]) {
                    edge_data$yend[j] <- edge_data$yend[j] + offset[3]
                    edge_data$y[paired_idx[j]] <- edge_data$y[paired_idx[j]] + offset[2]
                } else {
                    edge_data$yend[j] <- edge_data$yend[j] + offset[2]
                    edge_data$y[paired_idx[j]] <- edge_data$y[paired_idx[j]] + offset[3]
                }
            }
        }
    }

    edge_data$from <- NULL
    edge_data$to <- NULL
    edge_data
}


generate_label_data <- function(virus_info, hex_data) {
    x <- virus_info$x
    y <- virus_info$y

    width <- vapply(hex_data, function(d) max(d$x), numeric(1)) - x
    height <- vapply(hex_data, function(d) max(d$y), numeric(1)) - y

    xadj <- diff(range(virus_info$x)) * 0.02
    yadj <- diff(range(virus_info$y)) * 0.02
    hjust <- vjust <- 0.5

    i <- virus_info$label_position == "left"
    if (any(i)) {
        x[i] <- x[i] - width[i] - xadj
        hjust[i] <- 1
    }

    i <- virus_info$label_position == "right"
    if (any(i)) {
        x[i] <- x[i] + width[i] + xadj
        hjust[i] <- 0
    }

    i <- virus_info$label_position == "below"
    if (any(i)) {
        y[i] <- y[i] - height[i] - yadj
        vjust[i] <- 1
    }

    i <- virus_info$label_position == "upper"
    if (any(i)) {
        y[i] <- y[i] + height[i] + yadj
        vjust[i] <- 0
    }

    d <- data.frame(
        x = x,
        y = y,
        label = virus_info$label,
        vjust = vjust,
        hjust = hjust,
        stringsAsFactors = FALSE
    )
    d[!is.na(virus_info$label_position) & virus_info$label_position != "none", , drop = FALSE]
}


generate_elbow_flow_data <- function(flow_data, link_elbow_position = 0.5) {
    validate_link_elbow_position(link_elbow_position)

    xmid <- flow_data$x + (flow_data$xend - flow_data$x) * link_elbow_position

    lead <- data.frame(
        x = flow_data$x,
        xend = xmid,
        y = flow_data$y,
        yend = flow_data$y
    )
    middle <- data.frame(
        x = xmid,
        xend = xmid,
        y = flow_data$y,
        yend = flow_data$yend
    )
    tail <- data.frame(
        x = xmid,
        xend = flow_data$xend,
        y = flow_data$yend,
        yend = flow_data$yend
    )

    list(lead = lead, middle = middle, tail = tail)
}


estimate_asp <- function(ASP) {
    if (ASP < 1) {
        asp.x <- ASP
        asp.y <- 1
    } else {
        asp.x <- 1
        asp.y <- 1 / ASP
    }
    c(asp.x, asp.y)
}


asp_ <- function(virus_info, asp = 1) {
    diff(range(virus_info$x, na.rm = TRUE)) / diff(range(virus_info$y, na.rm = TRUE)) / asp
}
