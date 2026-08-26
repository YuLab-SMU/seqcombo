test_that("set_layout preserves x coordinates for character ids", {
    timeline_segments <- data.frame(
        id = rep(c("avian_1990", "human_1990", "swine_2000", "human_2009"), each = 3),
        sample_time = rep(c(1990, 1990, 2000, 2009), each = 3),
        segment = rep(c("PB2", "HA", "NA"), 4),
        color = c(
            rep("purple", 3),
            rep("red", 3),
            c("darkgreen", "red", "purple"),
            c("darkgreen", "lightgreen", "purple")
        ),
        host = rep(c("Avian", "Human", "Swine", "Human"), each = 3),
        stringsAsFactors = FALSE
    )

    flow_info <- build_flow_info(
        from = c("avian_1990", "human_1990", "swine_2000"),
        to = c("swine_2000", "swine_2000", "human_2009")
    )
    virus_info <- build_virus_info_from_long(
        timeline_segments,
        x = "sample_time",
        keep = "host"
    )

    layouted <- set_layout(virus_info, flow_info, preserve_x = TRUE)

    expect_equal(layouted$x, virus_info$x)
    expect_false(anyNA(layouted$y))
    expect_setequal(layouted$id, virus_info$id)
})

test_that("layout_timeline keeps time on the chosen axis", {
    timeline_segments <- data.frame(
        id = rep(c("avian_1990", "swine_2000", "human_2009"), each = 3),
        sample_time = rep(c(1990, 2000, 2009), each = 3),
        segment = rep(c("PB2", "HA", "NA"), 3),
        color = c(rep("purple", 3), rep("darkgreen", 3), rep("red", 3)),
        stringsAsFactors = FALSE
    )

    flow_info <- build_flow_info(
        from = c("avian_1990", "swine_2000"),
        to = c("swine_2000", "human_2009")
    )

    virus_x <- build_virus_info_from_long(timeline_segments, x = "sample_time")
    layout_x <- layout_timeline(virus_x, flow_info, time_col = "x", axis = "x")

    expect_equal(layout_x$x, virus_x$x)
    expect_false(anyNA(layout_x$y))

    virus_y <- build_virus_info_from_long(timeline_segments, y = "sample_time")
    layout_y <- layout_timeline(virus_y, flow_info, time_col = "y", axis = "y")

    expect_equal(layout_y$y, virus_y$y)
    expect_false(anyNA(layout_y$x))
})

test_that("set_layout rejects preserving both axes", {
    data <- example_seqcombo_data("basic")

    expect_error(
        set_layout(
            data$virus_info,
            data$flow_info,
            preserve_x = TRUE,
            preserve_y = TRUE
        ),
        "'preserve_x' and 'preserve_y' cannot both be TRUE"
    )
})
