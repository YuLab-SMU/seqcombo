##' Visualize virus reassortment events
##'
##' @title hybrid_plot
##' @param virus_info virus information
##' @param flow_info flow information
##' @param v_color the color of outer boundary of virus; can use expression
##'   (e.g. `v_color = ~Host`) to color virus by specific variable
##' @param v_fill the color to fill viruses; can use expression
##'   (e.g. `v_fill = ~Host`) to fill virus by specific variable
##' @param v_shape one of `hexagon` or `ellipse`
##' @param l_color color of the lines that indicate genetic flow
##' @param asp aspect ratio of the plotting device
##' @param parse whether parse label, only works if `label` and
##'   `label_position` exist
##' @param g_height height of regions to plot gene segments relative to the virus
##' @param g_width width of gene segment relative to width of the virus
##' @param t_size size of text label
##' @param t_color color of text label
##' @param link_style one of `segment`, `curve`, or `elbow`
##' @param link_curvature curvature used when `link_style = "curve"`
##' @param link_elbow_position horizontal bend position used when
##'   `link_style = "elbow"`
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom rlang .data
##' @export
##' @examples
##' n <- 8
##' virus_info <- build_virus_info(
##'     id = 1:7,
##'     x = c(rep(1990, 4), rep(2000, 2), 2009),
##'     y = c(1, 2, 3, 5, 1.5, 3, 4),
##'     segment_color = list(
##'         rep("purple", n),
##'         rep("red", n),
##'         rep("darkgreen", n),
##'         rep("lightgreen", n),
##'         c("darkgreen", "darkgreen", "red", "darkgreen", "red", "purple", "red", "purple"),
##'         c("darkgreen", "darkgreen", "red", "darkgreen", "darkgreen", "purple", "red", "purple"),
##'         c("darkgreen", "lightgreen", "lightgreen", "darkgreen",
##'           "darkgreen", "purple", "red", "purple")
##'     )
##' )
##' flow_info <- build_flow_info(
##'     from = c(1, 2, 3, 3, 4, 5, 6),
##'     to = c(5, 5, 5, 6, 7, 6, 7)
##' )
##' hybrid_plot(virus_info, flow_info, link_style = "elbow")
##' @author Guangchuang Yu
hybrid_plot <- function(virus_info, flow_info, v_color = "darkgreen",
                        v_fill = "steelblue", v_shape = "ellipse",
                        l_color = "black", asp = 1, parse = FALSE,
                        g_height = 0.65, g_width = 0.65, t_size = 3.88,
                        t_color = "black", link_style = "segment",
                        link_curvature = 0.15, link_elbow_position = 0.5) {
    ggplot(virus_info, aes(x = .data[["x"]], y = .data[["y"]])) +
        geom_hybrid(
            virus_info = virus_info,
            flow_info = flow_info,
            v_color = v_color,
            v_fill = v_fill,
            v_shape = v_shape,
            l_color = l_color,
            asp = asp,
            parse = parse,
            g_height = g_height,
            g_width = g_width,
            t_size = t_size,
            t_color = t_color,
            link_style = link_style,
            link_curvature = link_curvature,
            link_elbow_position = link_elbow_position
        )
}


##' Geom layer of genotype
##'
##' @title geom_genotype
##' @param virus_info virus information
##' @param v_color the color of outer boundary of virus; can use expression
##'   (e.g. `v_color = ~Host`) to color virus by specific variable
##' @param v_fill the color to fill viruses; can use expression
##'   (e.g. `v_fill = ~Host`) to fill virus by specific variable
##' @param v_shape one of `hexagon` or `ellipse`
##' @param l_color color of the lines that indicate genetic flow
##' @param asp aspect ratio of the plotting device
##' @param g_height height of regions to plot gene segments relative to the virus
##' @param g_width width of gene segment relative to width of the virus
##' @return geom layer
##' @export
##' @examples
##' n <- 8
##' virus_info <- build_virus_info(
##'     id = 1:7,
##'     x = c(rep(1990, 4), rep(2000, 2), 2009),
##'     y = c(1, 2, 3, 5, 1.5, 3, 4),
##'     segment_color = list(
##'         rep("purple", n),
##'         rep("red", n),
##'         rep("darkgreen", n),
##'         rep("lightgreen", n),
##'         c("darkgreen", "darkgreen", "red", "darkgreen", "red", "purple", "red", "purple"),
##'         c("darkgreen", "darkgreen", "red", "darkgreen", "darkgreen", "purple", "red", "purple"),
##'         c("darkgreen", "lightgreen", "lightgreen", "darkgreen",
##'           "darkgreen", "purple", "red", "purple")
##'     )
##' )
##' ggplot2::ggplot() + geom_genotype(virus_info)
##' @author Guangchuang Yu
geom_genotype <- function(virus_info, v_color = "darkgreen",
                          v_fill = "steelblue", v_shape = "ellipse",
                          l_color = "black", asp = 1, g_height = 0.65,
                          g_width = 0.65) {
    scene <- prepare_reassortment_scene(
        virus_info = virus_info,
        flow_info = NULL,
        asp = asp,
        v_shape = v_shape
    )

    build_genotype_layers(
        scene = scene,
        v_color = v_color,
        v_fill = v_fill,
        g_height = g_height,
        g_width = g_width
    )
}


##' Geom layer for reassortment events
##'
##' @title geom_hybrid
##' @inheritParams hybrid_plot
##' @return geom layer
##' @export
##' @examples
##' n <- 8
##' virus_info <- build_virus_info(
##'     id = 1:7,
##'     x = c(rep(1990, 4), rep(2000, 2), 2009),
##'     y = c(1, 2, 3, 5, 1.5, 3, 4),
##'     segment_color = list(
##'         rep("purple", n),
##'         rep("red", n),
##'         rep("darkgreen", n),
##'         rep("lightgreen", n),
##'         c("darkgreen", "darkgreen", "red", "darkgreen", "red", "purple", "red", "purple"),
##'         c("darkgreen", "darkgreen", "red", "darkgreen", "darkgreen", "purple", "red", "purple"),
##'         c("darkgreen", "lightgreen", "lightgreen", "darkgreen",
##'           "darkgreen", "purple", "red", "purple")
##'     )
##' )
##' flow_info <- build_flow_info(
##'     from = c(1, 2, 3, 3, 4, 5, 6),
##'     to = c(5, 5, 5, 6, 7, 6, 7)
##' )
##' ggplot2::ggplot() + geom_hybrid(virus_info, flow_info)
##' @author Guangchuang Yu
geom_hybrid <- function(virus_info, flow_info, v_color = "darkgreen",
                        v_fill = "steelblue", v_shape = "ellipse",
                        l_color = "black", asp = 1, parse = FALSE,
                        g_height = 0.65, g_width = 0.65, t_size = 3.88,
                        t_color = "black", link_style = "segment",
                        link_curvature = 0.15, link_elbow_position = 0.5) {
    scene <- prepare_reassortment_scene(
        virus_info = virus_info,
        flow_info = flow_info,
        asp = asp,
        v_shape = v_shape
    )

    c(
        build_genotype_layers(
            scene = scene,
            v_color = v_color,
            v_fill = v_fill,
            g_height = g_height,
            g_width = g_width
        ),
        list(
            build_flow_layer(
                scene = scene,
                l_color = l_color,
                link_style = link_style,
                link_curvature = link_curvature,
                link_elbow_position = link_elbow_position
            ),
            build_label_layer(
                scene = scene,
                parse = parse,
                t_size = t_size,
                t_color = t_color
            )
        )
    )
}


##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 geom_curve
##' @importFrom ggplot2 geom_text
##' @importFrom grid unit
##' @importFrom grid arrow
##' @importFrom yulab.utils get_fun_from_pkg
build_genotype_layers <- function(scene, v_color, v_fill, g_height, g_width) {
    default_aes <- aes(x = .data[["x"]], y = .data[["y"]])

    virus_capsule <- geom_virus_capsule(
        mapping = default_aes,
        virus_info = scene$virus_info,
        hex_data = scene$capsule_data,
        color = v_color,
        fill = v_fill
    )

    virus_segment <- lapply(seq_len(nrow(scene$virus_info)), function(i) {
        geom_gene_segment(
            hexd = scene$capsule_data[[i]],
            color = scene$virus_info$segment_color[[i]],
            g_height = g_height,
            g_width = g_width
        )
    })

    list(virus_capsule, virus_segment)
}


build_flow_layer <- function(scene, l_color, link_style, link_curvature,
                             link_elbow_position) {
    if (is.null(scene$flow_data)) {
        return(NULL)
    }

    link_style <- validate_link_style(link_style)
    mapping <- aes(
        x = .data[["x"]],
        xend = .data[["xend"]],
        y = .data[["y"]],
        yend = .data[["yend"]]
    )

    if (link_style == "curve") {
        return(geom_curve(
            mapping = mapping,
            data = scene$flow_data,
            arrow = arrow(length = unit(0.3, "cm")),
            color = l_color,
            curvature = link_curvature,
            inherit.aes = FALSE
        ))
    }

    if (link_style == "elbow") {
        elbow_data <- generate_elbow_flow_data(scene$flow_data, link_elbow_position)
        return(list(
            geom_segment(
                mapping = mapping,
                data = rbind(elbow_data$lead, elbow_data$middle),
                color = l_color,
                inherit.aes = FALSE
            ),
            geom_segment(
                mapping = mapping,
                data = elbow_data$tail,
                arrow = arrow(length = unit(0.3, "cm")),
                color = l_color,
                inherit.aes = FALSE
            )
        ))
    }

    geom_segment(
        mapping = mapping,
        data = scene$flow_data,
        arrow = arrow(length = unit(0.3, "cm")),
        color = l_color,
        inherit.aes = FALSE
    )
}


build_label_layer <- function(scene, parse, t_size, t_color) {
    if (is.null(scene$label_data) || nrow(scene$label_data) == 0) {
        return(NULL)
    }

    ld <- scene$label_data
    if (identical(parse, "emoji")) {
        emoji <- get_fun_from_pkg("emojifont", "emoji")
        ld$label <- emoji(ld$label)
        ld$vjust <- ld$vjust - 0.25
        parse <- FALSE
        family <- "EmojiOne"
    } else {
        family <- "sans"
    }

    geom_text(
        aes(
            x = .data[["x"]],
            y = .data[["y"]],
            label = .data[["label"]],
            vjust = .data[["vjust"]],
            hjust = .data[["hjust"]]
        ),
        data = ld,
        parse = parse,
        family = family,
        size = t_size,
        color = t_color,
        inherit.aes = FALSE
    )
}
