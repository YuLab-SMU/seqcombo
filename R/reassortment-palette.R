##' Host color presets for seqcombo figures
##'
##' A curated color palette for common influenza A host origins, plus helpers
##' to map host labels to colors and to wire the palette into ggplot2 scales.
##'
##' @param levels character vector of host labels, or \code{NULL} for the
##' full built-in set
##' @return a named character vector of colors
##' @examples
##' seqcombo_host_palette()
##' seqcombo_host_palette(c("avian", "human"))
##' @author Guangchuang Yu
##' @export
seqcombo_host_palette <- function(levels = NULL) {
    builtin <- .seqcombo_host_defaults()

    if (is.null(levels)) {
        return(builtin)
    }

    levels <- as.character(levels)
    lookup <- stats::setNames(
        seq_along(builtin),
        tolower(names(builtin))
    )
    hit <- tolower(levels) %in% names(lookup)

    out <- rep(NA_character_, length(levels))
    idx <- match(tolower(levels[hit]), names(lookup))
    out[hit] <- unname(builtin)[idx]

    extra_n <- sum(!hit)
    if (extra_n > 0) {
        extra <- grDevices::hcl(
            h = seq(15, 375, length.out = extra_n + 1)[-1],
            c = 78,
            l = 55
        )
        out[!hit] <- unname(extra)
    }

    stats::setNames(out, levels)
}


##' Apply the seqcombo host palette to a vector of labels
##'
##' Map host labels (or any categorical labels) to colors using the
##' [seqcombo_host_palette()] presets, which is handy for filling the
##' \code{color} column consumed by \code{segment_color}.
##'
##' @param x character vector of labels
##' @param levels optional character vector controlling level order and
##' color assignment; defaults to unique values of \code{x} in order of
##' appearance
##' @return an unnamed character vector of colors with one entry per element
##' in \code{x}
##' @examples
##' hosts <- c("Avian", "Avian", "Human", "Swine")
##' apply_seqcombo_palette(hosts)
##' @author Guangchuang Yu
##' @export
apply_seqcombo_palette <- function(x, levels = NULL) {
    x <- as.character(x)

    if (is.null(levels)) {
        levels <- unique(x[!is.na(x)])
    } else {
        levels <- as.character(levels)
    }

    pal <- seqcombo_host_palette(levels)
    cols <- unname(pal)[match(tolower(x), tolower(names(pal)))]
    cols[is.na(cols)] <- "#BEBEBE"
    cols
}


##' ggplot2 scales using the seqcombo host palette
##'
##' Create paired colour/fill discrete scales wired to
##' [seqcombo_host_palette()], so plots colored by host metadata need no
##' hand-written \code{scale_color_manual()} / \code{scale_fill_manual()}
##' calls.
##'
##' @param aesthetics character vector selecting which aesthetics to supply,
##' any of \code{"colour"} (\code{"color"}) and \code{"fill"}
##' @param ... additional arguments passed to
##' \code{\link[ggplot2]{scale_colour_manual}}
##' @return a list of ggplot2 scales that can be added to a plot
##' @examples
##' data <- example_seqcombo_data()
##' autoplot(data, v_color = ~Host, v_fill = ~Host) +
##'     scale_seqcombo_host()
##' @author Guangchuang Yu
##' @export
scale_seqcombo_host <- function(aesthetics = c("colour", "fill"), ...) {
    aesthetics <- match.arg(aesthetics,
                            choices = c("colour", "color", "fill"),
                            several.ok = TRUE)

    values <- seqcombo_host_palette()
    lapply(aesthetics, function(aes_name) {
        if (aes_name == "fill") {
            ggplot2::scale_fill_manual(values = values, ...)
        } else {
            ggplot2::scale_colour_manual(values = values, ...)
        }
    })
}


.seqcombo_host_defaults <- function() {
    c(
        Human = "#E64B35",
        Avian = "#4DBBD5",
        Swine = "#00A087",
        Equine = "#3C5488",
        Bat = "#91D1C2"
    )
}
