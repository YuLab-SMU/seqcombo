##' Create a seqcombo plotting object
##'
##' Bundle `virus_info` and `flow_info` into a lightweight `seqcombo_data`
##' object for validation and autoplot methods.
##'
##' @param virus_info virus information
##' @param flow_info optional flow information
##' @return an object of class `seqcombo_data`
##' @examples
##' data <- example_seqcombo_data()
##' seqcombo_data <- as_seqcombo_data(data$virus_info, data$flow_info)
##' class(seqcombo_data)
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
##' Return ready-to-plot `seqcombo_data` objects covering common teaching
##' scenarios.
##'
##' @param type one of `basic`, `timeline`, `long`, or `genotype`. See Details.
##' @return an object of class `seqcombo_data`
##'
##' @details Available types:
##' \describe{
##'   \item{basic}{the classic hybrid reassortment network (numeric ids)}
##'   \item{timeline}{same biology arranged by `layout_timeline()`}
##'   \item{long}{built from tidy segment/flow tables via
##'     `build_virus_info_from_long()` and `build_flow_info_from_long()`; raw
##'     tables are attached as `attr(data, "long_tables")`}
##'   \item{genotype}{two parents plus two reassortants, without flows, for
##'     genotype-only plots via `geom_genotype()`}
##' }
##' Segment colors in `long` and `genotype` encode ancestral host through
##' [apply_seqcombo_palette()], so they pair naturally with
##' [scale_seqcombo_host()].
##'
##' @examples
##' data <- example_seqcombo_data("long")
##' table(attr(data, "long_tables")$flows)
##' @export
##' @author Guangchuang Yu
example_seqcombo_data <- function(type = c("basic", "timeline", "long", "genotype")) {
    type <- match.arg(type)

    if (type == "basic") {
        return(.example_basic(with_layout = FALSE))
    }
    if (type == "timeline") {
        return(.example_basic(with_layout = TRUE))
    }
    if (type == "long") {
        return(.example_from_long())
    }

    .example_genotype()
}


.example_basic <- function(with_layout = FALSE) {
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

    if (with_layout) {
        virus_info <- layout_timeline(virus_info, flow_info, time_col = "x")
    }

    as_seqcombo_data(virus_info, flow_info)
}


.example_from_long <- function() {
    segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
    ids <- c("avian_1996", "swine_eurasian_2002",
             "swine_triple_2008", "pandemic_h1n1_2009")
    times <- c(1996, 2002, 2008, 2009)
    sampling_host <- c("Avian", "Swine", "Swine", "Human")

    ## ancestral host of every segment: color carries ancestry information
    segment_origin <- rbind(
        c("Avian", "Avian", "Avian", "Avian", "Avian", "Avian", "Avian", "Avian"),
        c("Avian", "Avian", "Avian", "Avian", "Avian", "Human", "Avian", "Avian"),
        c("Human", "Swine", "Avian", "Swine", "Avian", "Avian", "Swine", "Swine"),
        c("Swine", "Human", "Swine", "Swine", "Avian", "Swine", "Swine", "Swine")
    )

    segment_df <- data.frame(
        id = rep(ids, each = length(segments)),
        sample_time = rep(times, each = length(segments)),
        segment = segments,
        origin_host = as.vector(segment_origin),
        sampling_host = rep(sampling_host, each = length(segments)),
        stringsAsFactors = FALSE
    )
    segment_df$color <- apply_seqcombo_palette(segment_df$origin_host)

    flow_df <- data.frame(
        from = c("avian_1996", "avian_1996",
                 "swine_eurasian_2002", "swine_triple_2008"),
        to = c("swine_eurasian_2002", "swine_triple_2008",
               "swine_triple_2008", "pandemic_h1n1_2009"),
        support = c(0.62, 0.84, 0.58, 0.91),
        stringsAsFactors = FALSE
    )

    virus_info <- build_virus_info_from_long(
        segment_df,
        id = "id", segment = "segment", color = "color",
        x = "sample_time", keep = "sampling_host"
    )
    virus_info$label <- c("Avian reservoir", "Eurasian swine H1N1",
                          "Triple reassortant swine", "Pandemic H1N1")
    virus_info$label_position <- "left"
    flow_info <- build_flow_info_from_long(flow_df, weight = "support")
    virus_info <- layout_timeline(virus_info, flow_info, time_col = "x")

    data <- as_seqcombo_data(virus_info, flow_info)
    attr(data, "long_tables") <- list(segments = segment_df, flows = flow_df)
    data
}


.example_genotype <- function() {
    segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
    pal <- .seqcombo_host_defaults()
    avian <- unname(pal[["Avian"]])
    human <- unname(pal[["Human"]])

    virus_info <- build_virus_info(
        id = c("avian_parent", "human_parent", "reassortant_i", "reassortant_ii"),
        x = c(1, 5, 1.7, 4.3),
        y = c(1, 1, 2.8, 2.8),
        segment_color = list(
            rep(avian, length(segments)),
            rep(human, length(segments)),
            c(rep(human, 4), rep(avian, 4)),
            c(human, avian, avian, human, avian, human, human, human)
        ),
        segment_name = rep(list(segments), 4),
        label = c("Avian parent", "Human parent",
                  "Reassortant I", "Reassortant II"),
        label_position = "below"
    )

    as_seqcombo_data(virus_info)
}


##' Autoplot a seqcombo plotting object
##'
##' Create a reassortment plot directly from a `seqcombo_data` object.
##'
##' @param object a `seqcombo_data` object
##' @param ... additional parameters passed to `hybrid_plot()` or `geom_genotype()`
##' @return a ggplot object
##' @examples
##' data <- example_seqcombo_data("timeline")
##' ggplot2::autoplot(data, link_style = "curve")
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
