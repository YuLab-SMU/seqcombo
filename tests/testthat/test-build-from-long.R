test_that("build_virus_info_from_long aggregates segment records by virus", {
    segment_df <- data.frame(
        id = rep(c("avian", "human"), each = 3),
        segment = rep(c("HA", "PB2", "NA"), 2),
        color = c("purple", "purple", "purple", "red", "red", "red"),
        sample_time = rep(c(1990, 2009), each = 3),
        host = rep(c("Avian", "Human"), each = 3),
        stringsAsFactors = FALSE
    )

    virus_info <- build_virus_info_from_long(
        segment_df,
        x = "sample_time",
        keep = "host",
        segment_order = c("PB2", "HA", "NA")
    )

    expect_equal(virus_info$id, c("avian", "human"))
    expect_equal(virus_info$x, c(1990, 2009))
    expect_equal(virus_info$host, c("Avian", "Human"))
    expect_equal(virus_info$segment_name[[1]], c("PB2", "HA", "NA"))
    expect_equal(virus_info$segment_color[[2]], c("red", "red", "red"))
})

test_that("build_virus_info_from_long rejects incomplete segment order", {
    segment_df <- data.frame(
        id = rep("avian", 2),
        segment = c("HA", "PB2"),
        color = c("purple", "purple"),
        stringsAsFactors = FALSE
    )

    expect_error(
        build_virus_info_from_long(
            segment_df,
            segment_order = "HA"
        ),
        "all segment values must be present in 'segment_order'"
    )
})

test_that("build_flow_info_from_long aggregates edges and weights", {
    flow_long <- data.frame(
        from = c("avian", "avian", "human"),
        to = c("swine", "swine", "swine"),
        segment = c("HA", "NA", "PB2"),
        weight = c(1, 2, 4),
        lineage = c("alpha", "alpha", "beta"),
        stringsAsFactors = FALSE
    )

    flow_info <- build_flow_info_from_long(
        flow_long,
        segment = "segment",
        weight = "weight",
        keep = "lineage"
    )

    avian_to_swine <- flow_info[flow_info$from == "avian" & flow_info$to == "swine", ]
    human_to_swine <- flow_info[flow_info$from == "human" & flow_info$to == "swine", ]

    expect_equal(nrow(flow_info), 2)
    expect_equal(avian_to_swine$weight, 3)
    expect_equal(avian_to_swine$lineage, "alpha")
    expect_equal(avian_to_swine$segment_name[[1]], c("HA", "NA"))
    expect_equal(human_to_swine$weight, 4)
})
