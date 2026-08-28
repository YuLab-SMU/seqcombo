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

test_that("set_layout spreads viruses sharing the same preserved coordinate", {
    timeline_segments <- data.frame(
        id = rep(c("a_1990", "b_1990", "c_1990", "human_2009"), each = 3),
        sample_time = rep(c(1990, 1990, 1990, 2009), each = 3),
        segment = rep(c("PB2", "HA", "NA"), 4),
        color = c(
            rep("purple", 3),
            rep("red", 3),
            rep("darkgreen", 3),
            rep("blue", 3)
        ),
        stringsAsFactors = FALSE
    )

    flow_info <- build_flow_info(
        from = c("a_1990", "b_1990", "c_1990"),
        to = c("human_2009", "human_2009", "human_2009")
    )
    virus_info <- build_virus_info_from_long(timeline_segments, x = "sample_time")

    layouted <- set_layout(virus_info, flow_info, preserve_x = TRUE, spread = TRUE)

    expect_equal(layouted$x, virus_info$x)
    expect_equal(length(unique(layouted$y)), nrow(layouted))
    ## slots are integers with unit spacing within each time group
    expect_true(all(layouted$y == round(layouted$y)))
    ## later time groups sit below earlier ones
    expect_equal(layouted$y[layouted$id == "human_2009"], 1)
    expect_equal(sort(layouted$y[layouted$id != "human_2009"]), 3:5)
})

test_that("spread_grouped_values assigns slotted positions", {
    res <- seqcombo:::spread_grouped_values(c(10, 20, 30, 40), c("a", "a", "b", "c"))
    expect_equal(res, c(6, 5, 3, 1))

    ## missing groups are kept as their own block (placed on top)
    res2 <- seqcombo:::spread_grouped_values(c(1, 2), c(NA, "x"))
    expect_equal(res2, c(3, 1))
})

test_that("layout_timeline spreads shared time points by default", {
    timeline_segments <- data.frame(
        id = rep(c("a_1990", "b_1990", "human_2009"), each = 3),
        sample_time = rep(c(1990, 1990, 2009), each = 3),
        segment = rep(c("PB2", "HA", "NA"), 3),
        color = c(rep("purple", 3), rep("red", 3), rep("darkgreen", 3)),
        stringsAsFactors = FALSE
    )

    flow_info <- build_flow_info(
        from = c("a_1990", "b_1990"),
        to = c("human_2009", "human_2009")
    )
    virus_info <- build_virus_info_from_long(timeline_segments, x = "sample_time")

    layouted <- layout_timeline(virus_info, flow_info, time_col = "x", axis = "x")

    expect_equal(layouted$x, virus_info$x)
    expect_true(all(layouted$y == round(layouted$y)))
    expect_equal(length(unique(layouted$y)), nrow(layouted))
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
