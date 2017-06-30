##' visualize virus reassortment events
##'
##'
##' @title plot_reassort
##' @param virus_info virus information
##' @param flow_info flow information
##' @param v_color the color of outer boundary of virus; can use expression (e.g. v_color=~Host) to color virus by specific variable
##' @param v_fill the color to fill viruses; can use expression (e.g. v_fill=~Host) to fill virus by specific variable
##' @param l_color color of the lines that indicate genetic flow
##' @param asp aspect ratio of the plotting device
##' @param parse whether parse label, only works if 'label' and 'label_position' exist
##' @param t_size size of text label
##' @param t_color color of text label
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 geom_blank
##' @importFrom ggplot2 aes_
##' @importFrom grid unit
##' @importFrom grid arrow
##' @importFrom rvcheck get_fun_from_pkg
##' @export
##' @examples
##' library(tibble)
##' n <- 8
##' virus_info <- tibble(id = 1:7,
##' x = c(rep(1990, 4), rep(2000, 2), 2009),
##' y = c(1,2,3,5, 1.5, 3, 4),
##' segment_color = list(rep('purple', n),
##' rep('red', n), rep('darkgreen', n), rep('lightgreen', n),
##' c('darkgreen', 'darkgreen', 'red', 'darkgreen', 'red', 'purple', 'red', 'purple'),
##' c('darkgreen', 'darkgreen', 'red', 'darkgreen', 'darkgreen', 'purple', 'red', 'purple'),
##' c('darkgreen', 'lightgreen', 'lightgreen', 'darkgreen', 'darkgreen', 'purple', 'red', 'purple')))
##'
##' flow_info <- tibble(from = c(1,2,3,3,4,5,6), to = c(5,5,5,6,7,6,7))
##'
##' plot_reassort(virus_info, flow_info)
##'
##' @author guangchuang yu
plot_reassort <- function(virus_info, flow_info, v_color="darkgreen", v_fill="steelblue",
                          l_color="black", asp=1, parse=FALSE, t_size=3.88, t_color="black") {

    require_col <- c('x', 'y', 'id', 'segment_color')
    if (!all(require_col %in% colnames(virus_info)))
        stop("'x', 'y', and 'id' columns are required in 'virus_info'...")


    if (!all(c('from', 'to') %in% colnames(flow_info)))
        stop("'from' and 'to' columns are required in 'flow_info'...")

    if (!'virus_size' %in% colnames(virus_info))
        virus_info$virus_size <- 1

    ASP <-  diff(range(virus_info$x)) / diff(range(virus_info$y)) / asp
    hex_data <- lapply(1:nrow(virus_info), function(i) {
        x <- generate_hex_data(
            x = virus_info$x[i],
            y = virus_info$y[i],
            size = virus_info$virus_size[i],
            ASP = ASP)
        x$id <- virus_info$id[i]
        return(x)
    })

    default_aes <- aes_(x=~x, y=~y)

    virus_capsule <- geom_virus_capsule(default_aes, virus_info, hex_data, v_color, v_fill, ASP)

    virus_segment <- lapply(1:nrow(virus_info), function(i)
        geom_gene_segment(hex_data[[i]],
                  color=virus_info$segment_color[[i]],
                  width = .65)
        )


    d <- generate_segment_data(virus_info, flow_info, hex_data, ASP)
    virus_link <- geom_segment(aes_(x=~x, xend=~xend, y=~y, yend=~yend), data=d, arrow=arrow(length=unit(.3, 'cm')), color=l_color)


    p <- ggplot(virus_info, default_aes) + geom_blank() + virus_capsule + virus_segment + virus_link

    if (all(c('label', 'label_position') %in% colnames(virus_info))) {
        ld <- generate_label_data(virus_info, hex_data)

        if (parse == 'emoji') {
            emoji <- get_fun_from_pkg("emojifont", "emoji")
            ld$label <- emoji(ld$label)
            ld$vjust <- ld$vjust - 0.25
            parse <- FALSE
            family <- "EmojiOne"
        } else {
            family <- 'sans'
        }

        p <- p + geom_text(aes_(x=~x, y=~y, label=~label, vjust=~vjust, hjust=~hjust),
                       data=ld, parse=parse, family=family, size=t_size, color=t_color, inherit.aes=FALSE)
    }

    return(p)
}

##' @importFrom ggplot2 geom_polygon
##' @importFrom utils modifyList
geom_virus_capsule <- function(mapping, virus_info, hex_data, color, fill, ASP=1, alpha=0.5, size=1) {
    hex.df <- do.call('rbind', hex_data)

    if (typeof(color) == "language") {
        vcol <- all.vars(color)
        if (!vcol %in% colnames(virus_info)) {
            stop("color variable not available...")
        }
        hex.df[, vcol] <- virus_info[[vcol]][match(hex.df$id, virus_info$id)]
        mapping <- modifyList(mapping, aes_(color=color))
    }

    if (typeof(fill) == "language") {
        vf <- all.vars(fill)
        if (!vf %in% colnames(virus_info)) {
            stop("fill variable not available...")
        }
        hex.df[,vf] <- virus_info[[vf]][match(hex.df$id, virus_info$id)]
        mapping <- modifyList(mapping, aes_(fill=fill))
    }
    mapping <- modifyList(mapping, aes_(group=~id))

    params <- list(mapping=mapping, data=hex.df, alpha=alpha, size=size, inherit.aes=FALSE)
    if (typeof(color) == "character") {
        params <- modifyList(list(color=color), params)
    }
    if (typeof(fill) == "character") {
        params <- modifyList(list(fill=fill), params)
    }

    do.call(geom_polygon, params)
}

##' @importFrom ggplot2 geom_rect
geom_gene_segment <- function(hexd, color, width=0.68) {
    n <- length(color)
    y <- hexd$y
    y <- y[y != min(y) & y != max(y)]
    ymin <- ymax <- seq(min(y), max(y), length.out=n+1)
    ymin <- ymin[-(n+1)]
    ymax <- ymax[-1]
    adjust <- (ymax - ymin) * (1 - width)/2
    ymin <- ymin + adjust
    ymax <- ymax - adjust

    xmin <- min(hexd$x)
    xmax <- max(hexd$x)
    xadj <- (xmax - xmin) * .1

    d <- data.frame(xmin = xmin + xadj,
                    xmax = xmax - xadj,
                    ymin = ymin,
                    ymax = ymax,
                    color = color)

    ## ## cannot coexists with aes(fill=VAR)
    ## geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=I(color)),
    ##           data= d, inherit.aes=F, show.legend=FALSE)

    dd <- split(d, color)
    lapply(seq_along(dd), function(i)
        geom_rect(aes_(xmin=~xmin, ymin=~ymin, xmax=~xmax, ymax=~ymax),
                  data= dd[[i]], fill=dd[[i]]$color[1], inherit.aes=F, show.legend=FALSE)
        )
}


generate_hex_data <- function(x, y, size, ASP) {
    asp <- estimate_asp(ASP)
    data.frame(x = c(rep(-sqrt(3)/2, 2), 0, rep(sqrt(3)/2, 2), 0) * size * asp[1] + x,
               y = c(0.5, -0.5, -1, -0.5, 0.5, 1)* size * asp[2] + y)
}


generate_segment_data <- function(virus_info, flow_info, hex_data, ASP=1) {
    x <- virus_info$x
    y <- virus_info$y

    x.from <- x[match(flow_info$from, virus_info$id)]
    x.to <- x[match(flow_info$to, virus_info$id)]
    y.from <- y[match(flow_info$from, virus_info$id)]
    y.to <- y[match(flow_info$to, virus_info$id)]

    width <- sapply(hex_data, function(x) max(x$x)) - x
    height <- sapply(hex_data, function(x) max(x$y)) - y
    names(width) <- names(height) <- virus_info$id

    xadj <- diff(range(virus_info$x)) * 0.01
    yadj <- diff(range(virus_info$y)) * 0.01

    asp <- estimate_asp(ASP)
    xdiff <- (x.to - x.from) * asp[1]
    ydiff <- (y.to - y.from) * asp[2]
    idx <- abs(ydiff) - abs(xdiff) > 0

    x.from.adj <- width[flow_info$from] + xadj
    x.from.adj[idx] <- 0
    x.to.adj <- width[flow_info$to] + xadj
    x.to.adj[idx] <- 0

    y.from.adj <- rep(0, nrow(flow_info))
    y.from.adj[idx] <- height[flow_info$from][idx] + yadj
    y.to.adj <- rep(0, nrow(flow_info))
    y.to.adj[idx] <- height[flow_info$to][idx] + yadj

    ii <- x.from != x.to
    x.direction <- sign(x.to - x.from)
    x.from <- x.from + x.from.adj * x.direction
    x.to <- x.to - x.to.adj * x.direction

    y.direction <- sign(y.to - y.from)
    y.from <- y.from + y.from.adj * y.direction
    y.to <- y.to - y.to.adj  * y.direction

    d <- data.frame(x = x.from, xend=x.to,
                    y = y.from, yend=y.to)

    xyadjust <- function(d, x, y, width, height, id) {
        d$id <- id
        ii <- which(duplicated(d[, c(x, y, 'id')]))
        if (length(ii)) {
            ii <- ii[!duplicated(d[ii,'id'])]

            for (i in ii) {
                j <- d$id == d$id[i] & d[,x] == d[i, x] & d[,y] == d[i, y]
                if (all(d$xend[j] > d$x[j])) {
                    ## left to right, adjust y
                    h <- height[d$id[j]][1]/2
                    offset <- seq(-h, h, length.out=sum(j)+2)
                    offset <- offset[-c(1, length(offset))]
                    if (y == 'y') {
                        jj <- order(d[j, 'yend'], decreasing=FALSE)
                    } else {
                        jj <- order(d[j, 'y'], decreasing=FALSE)
                    }
                    d[j,y] <- d[j,y] + offset[jj]
                } else {
                    ## adjust x
                    w <- width[d$id[j]][1]/2
                    offset <- seq(-w, w, length.out=sum(j)+2)
                    offset <- offset[-c(1, length(offset))]
                    if (y == 'x') {
                        jj <- order(d[j, 'xend'], decreasing=FALSE)
                    } else {
                        jj <- order(d[j, 'x'], decreasing=FALSE)
                    }
                    d[j,x] <- d[j,x] + offset[jj]
                }
            }
        }
        return(d)
    }

    d <- xyadjust(d, 'x', 'y', width, height, flow_info$from)
    d <- xyadjust(d, 'xend', 'yend', width, height, flow_info$to)
    return(d)
}


generate_label_data <- function(virus_info, hex_data) {
    x <- virus_info$x
    y <- virus_info$y

    width <- sapply(hex_data, function(x) max(x$x)) - x
    height <- sapply(hex_data, function(x) max(x$y)) - y

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

    d <- data.frame(x=x, y=y, label=virus_info$label, vjust=vjust, hjust=hjust)
    d <- d[virus_info$label_position != 'none',]

    return(d)
}



estimate_asp <- function(ASP) {
    if (ASP < 1) {
        asp.x <- ASP
        asp.y <- 1
    } else {
        asp.x <- 1
        asp.y <- 1/ASP
    }
    return(c(asp.x, asp.y))
}

