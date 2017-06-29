
hex_virus <- function(hexd, size, fill, color, segment_color,  width=.68, ...) {
    list(geom_polygon(aes(x, y), data=hexd, fill = fill, colour = color, size=size, inherit.aes=F, ...),
         geom_gene_segment(hexd, segment_color, width)
         )
}



geom_gene_segment <- function(hexd, color=rainbow(8), width=0.68) {
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
    xadj <- (xmax - xmin) * .075

    d <- data.frame(xmin = xmin + xadj,
                    xmax = xmax - xadj,
                    ymin = ymin,
                    ymax = ymax,
                    color = color)

    geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=I(color)), d, inherit.aes=F, show.legend=F)
}


generate_hex_data <- function(x, y, size, ASP) {
    data.frame(x = c(rep(-sqrt(3)/2, 2), 0, rep(sqrt(3)/2, 2), 0) * size + x,
               y = c(0.5, -0.5, -1, -0.5, 0.5, 1)/ASP * size + y)
}


generate_segment_data <- function(virus_info, flow_info, hex_data) {
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

    ii <- x.from != x.to

    x.direction <- sign(x.to - x.from)[ii]
    x.from[ii] <- x.from[ii] + (width[flow_info$from][ii] + xadj) * x.direction
    x.to[ii] <- x.to[ii] - (width[flow_info$to][ii] + xadj) * x.direction


    y.direction <- sign(y.to - y.from)[!ii]
    y.from[!ii] <- y.from[!ii] + (height[flow_info$from][!ii] + yadj) * y.direction
    y.to[!ii] <- y.to[!ii] - (height[flow_info$to][!ii] + yadj) * y.direction

    data.frame(x = x.from, xend=x.to,
               y = y.from, yend=y.to)

}


generate_label_data <- function(virus_info, hex_data) {
    x <- virus_info$x
    y <- virus_info$y

    width <- sapply(hex_data, function(x) max(x$x)) - x
    height <- sapply(hex_data, function(x) max(x$y)) - y

    xadj <- diff(range(virus_info$x)) * 0.05
    yadj <- diff(range(virus_info$y)) * 0.05


    i <- virus_info$label_position == "left"
    if (any(i))
        x[i] <- x[i] - width[i] - xadj

    i <- virus_info$label_position == "right"
    if (any(i))
        x[i] <- x[i] + width[i] + xadj

    i <- virus_info$label_position == "below"
    if (any(i))
        y[i] <- y[i] - height[i] - yadj

    i <- virus_info$label_position == "upper"
    if (any(i))
        y[i] <- y[i] + height[i] + yadj

    d <- data.frame(x=x, y=y, label=virus_info$label)
    d <- d[virus_info$label_position != 'none',]

    return(d)
}


plot_reassort <- function(virus_info, flow_info, fill="steelblue") {
    ASP <- diff(range(virus_info$x))/diff(range(virus_info$y))

    hex_data <- lapply(1:nrow(virus_info), function(i)
        generate_hex_data(
            x = virus_info$x[i],
            y = virus_info$y[i],
            size = virus_info$virus_size[i],
            ASP = ASP))

    viruses <- lapply(1:nrow(virus_info), function(i)
        hex_virus(hex_data[[i]],
                  fill=fill, color=virus_info$color[i],
                  segment_color=virus_info$segment_color[[i]],
                  size=1)
        )
    d <- generate_segment_data(virus_info, flow_info, hex_data)

    p <- ggplot(virus_info, aes(x, y)) + geom_blank() + viruses +
        geom_segment(aes(x=x, xend=xend, y=y, yend=yend), data=d, arrow=arrow(length=unit(.3, 'cm')))

    if (all(c('label', 'label_position') %in% colnames(virus_info))) {
        ld <- generate_label_data(virus_info, hex_data)
        p <- p + geom_text(aes(x, y, label=label),
                           data=ld, inherit.aes=FALSE)
    }

    return(p)
}




require(tibble)
require(ggplot2)

n <- 8
virus_info <- tibble(id = 1:7,
                     x = c(rep(1990, 4), rep(2000, 2), 2009),
                     y = c(1,2,3,5, 1.5, 3, 3),
                     segment_color = list(rep('purple', n),
                                          rep('red', n),
                                          rep('darkgreen', n),
                                          rep('lightgreen', n),
                                          c('darkgreen', 'darkgreen', 'red', 'darkgreen', 'red', 'purple', 'red', 'purple'),
                                          c('darkgreen', 'darkgreen', 'red', 'darkgreen', 'darkgreen', 'purple', 'red', 'purple'),
                                          c('darkgreen', 'lightgreen', 'lightgreen', 'darkgreen', 'darkgreen', 'purple', 'red', 'purple')),
                     color = c("purple", "red", rep("darkgreen", 3), "red", "darkgreen"),
                     virus_size = c(rep(1, 3), 2, 1, 1, 1.5),
                     label = c("Avian", "Human\nH3N2", "Classic\nswine\nH1N1", "Eurasian swine", rep("North American swin\n H3N2 and H1N2", 2), "2009 Human H1N1"),
                     label_position = c('left', 'left', 'left', 'below', 'below', 'none', 'below')
                     )

flow_info <- tibble(from = c(1,2,3,3,4,5,6),
                    to = c(5,5,5,6,7,6,7))

plot_reassort(virus_info, flow_info)


