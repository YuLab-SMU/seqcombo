##' Create a seqcombo plotting object
##'
##' Bundle `virus_info` and `flow_info` into a lightweight `seqcombo_data`
##' object for validation and autoplot methods.
##'
##' @param virus_info virus information
##' @param flow_info optional flow information
##' @return an object of class `seqcombo_data`
##' @export
##' @author Guangchuang Yu
as_seqcombo_data <- function(virus_info, flow_info = NULL) {
    structure(
        list(
            virus_info = virus_info,
            flow_info = flow_info
        ),
        class = "seqcombo_data"
    )
}


##' Check seqcombo plotting inputs
##'
##' Validate `virus_info` and `flow_info` before plotting.
##'
##' @param x a `seqcombo_data` object or a `virus_info` data frame
##' @param flow_info optional `flow_info` when `x` is a `virus_info` data frame
##' @param require_coordinates whether `x` and `y` coordinates must already exist
##' @return invisibly returns `TRUE` when the inputs are valid
##' @export
##' @author Guangchuang Yu
check_seqcombo_data <- function(x, flow_info = NULL, require_coordinates = TRUE) {
    data <- normalize_seqcombo_input(x, flow_info)

    validate_virus_info(data$virus_info, require_coordinates = require_coordinates)
    validate_segment_color(data$virus_info)
    validate_segment_name(data$virus_info)
    if (!is.null(data$flow_info)) {
        validate_flow_info(data$flow_info, data$virus_info$id)
    }

    invisible(TRUE)
}


##' Example seqcombo plotting data
##'
##' Return a ready-to-plot `seqcombo_data` object for demos and tests.
##'
##' @param type one of `basic` or `timeline`
##' @return an object of class `seqcombo_data`
##' @export
##' @author Guangchuang Yu
example_seqcombo_data <- function(type = c("basic", "timeline")) {
    type <- match.arg(type)

    n <- 8
    segment_name <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
    virus_info <- build_virus_info(
        id = 1:7,
        x = c(rep(1990, 4), rep(2000, 2), 2009),
        y = c(1, 2, 3, 5, 1.5, 3, 4),
        segment_color = list(
            rep("purple", n),
            rep("red", n),
            rep("darkgreen", n),
            rep("lightgreen", n),
            c("darkgreen", "darkgreen", "red", "darkgreen", "red", "purple", "red", "purple"),
            c("darkgreen", "darkgreen", "red", "darkgreen", "darkgreen", "purple", "red", "purple"),
            c("darkgreen", "lightgreen", "lightgreen", "darkgreen",
              "darkgreen", "purple", "red", "purple")
        ),
        segment_name = rep(list(segment_name), 7),
        Host = c("Avian", "Human", rep("Swine", 4), "Human"),
        label = c("Avian", "Human\nH3N2", "Classic\nswine\nH1N1", "Eurasian swine",
                  "North American swine\ntriple reassortant H3N2",
                  "North American swine\ntriple reassortant H1N2", "2009 Human H1N1"),
        label_position = c("left", "left", "left", "below", "below", "upper", "below"),
        virus_size = c(rep(1, 3), 2, 1, 1, 1.5)
    )

    flow_info <- build_flow_info(
        from = c(1, 2, 3, 3, 4, 5, 6),
        to = c(5, 5, 5, 6, 7, 6, 7)
    )

    if (type == "timeline") {
        virus_info <- layout_timeline(virus_info, flow_info, time_col = "x")
    }

    as_seqcombo_data(virus_info, flow_info)
}


##' Autoplot a seqcombo plotting object
##'
##' Create a reassortment plot directly from a `seqcombo_data` object.
##'
##' @param object a `seqcombo_data` object
##' @param ... additional parameters passed to `hybrid_plot()` or `geom_genotype()`
##' @return a ggplot object
##' @importFrom ggplot2 autoplot
##' @export
##' @author Guangchuang Yu
autoplot.seqcombo_data <- function(object, ...) {
    check_seqcombo_data(object, require_coordinates = TRUE)

    if (is.null(object$flow_info)) {
        return(ggplot2::ggplot() + geom_genotype(object$virus_info, ...))
    }

    hybrid_plot(object$virus_info, object$flow_info, ...)
}


##' Publication-oriented theme for seqcombo figures
##'
##' A compact ggplot2 theme tuned for reassortment diagrams.
##'
##' @param base_size base text size
##' @param base_family base font family
##' @return a ggplot2 theme object
##' @export
##' @author Guangchuang Yu
theme_seqcombo_pub <- function(base_size = 11, base_family = "sans") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            legend.title = ggplot2::element_text(face = "bold"),
            legend.key.height = grid::unit(0.8, "lines"),
            plot.title = ggplot2::element_text(face = "bold"),
            plot.caption = ggplot2::element_text(hjust = 0)
        )
}


normalize_seqcombo_input <- function(x, flow_info = NULL) {
    if (inherits(x, "seqcombo_data")) {
        return(x)
    }

    as_seqcombo_data(x, flow_info)
}
