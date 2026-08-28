test_that("generate_elbow_flow_data carries metadata columns", {
    flow_data <- data.frame(
        x = c(0, 0), xend = c(10, 10),
        y = c(0, 5), yend = c(5, 0),
        alpha = c(0.2, 0.8), linewidth = c(1, 3),
        stringsAsFactors = FALSE
    )
    segs <- seqcombo:::generate_elbow_flow_data(flow_data, link_elbow_position = 0.5)

    for (part in c("lead", "middle", "tail")) {
        expect_true(all(c("alpha", "linewidth") %in% colnames(segs[[part]])))
        expect_equal(nrow(segs[[part]]), 2)
    }
    expect_equal(segs$lead$alpha, c(0.2, 0.8))
    expect_equal(segs$middle$y, c(0, 5))
    expect_equal(segs$middle$yend, c(5, 0))
    expect_equal(segs$tail$xend, c(10, 10))
})

test_that("generate_elbow_flow_data staggers shared elbow positions", {
    flow_data <- data.frame(
        x = rep(0, 4), xend = rep(10, 4),
        y = 1:4, yend = 1:4
    )
    segs <- seqcombo:::generate_elbow_flow_data(
        flow_data,
        link_elbow_position = 0.5,
        elbow_spread = 0.08
    )
    xmid <- segs$middle$x
    expect_equal(length(unique(xmid)), 4)
    expect_true(all(xmid > 0 & xmid < 10))
    expect_equal(mean(xmid), 5)
    expect_equal(segs$lead$xend, xmid)
    expect_equal(segs$tail$x, xmid)

    ## explicit per-flow positions are honoured as-is
    segs2 <- seqcombo:::generate_elbow_flow_data(
        flow_data,
        link_elbow_position = c(0.2, 0.4, 0.6, 0.8)
    )
    expect_equal(segs2$middle$x, c(2, 4, 6, 8))

    ## staggered positions are clamped inside (0, 1)
    segs3 <- seqcombo:::generate_elbow_flow_data(
        flow_data,
        link_elbow_position = 0.05,
        elbow_spread = 0.05
    )
    expect_lte(segs3$middle$x[1], 1e-6)
    segs4 <- seqcombo:::generate_elbow_flow_data(
        flow_data,
        link_elbow_position = 0.95,
        elbow_spread = 0.05
    )
    expect_gte(segs4$middle$x[4], 10 - 1e-6)

    expect_error(
        seqcombo:::generate_elbow_flow_data(flow_data, link_elbow_position = 1.5),
        "'link_elbow_position' must be numeric values between 0 and 1"
    )
    expect_error(
        seqcombo:::generate_elbow_flow_data(flow_data, link_elbow_position = c(0.5, 0.5)),
        "length of 'link_elbow_position' must be 1 or the number of flows"
    )
})

test_that("generate_label_data carries facet metadata and drops list columns", {
    virus_info <- data.frame(
        id = c("A", "B", "C"),
        x = c(0, 10, 20), y = c(0, 0, 0),
        label = c("virus A", "virus B", "virus C"),
        label_position = c("upper", "upper", "none"),
        Lineage = c("Europe", "NorthAmerica", "Asia"),
        segment_color = I(list(c("red"), c("blue"), c("green"))),
        stringsAsFactors = FALSE
    )
    hex_data <- list(
        A = data.frame(x = c(-1, 1), y = c(-1, 1)),
        B = data.frame(x = c(9, 11), y = c(-1, 1)),
        C = data.frame(x = c(19, 21), y = c(-1, 1))
    )

    ld <- seqcombo:::generate_label_data(virus_info, hex_data)

    expect_equal(ld$label, c("virus A", "virus B"))
    expect_true("Lineage" %in% colnames(ld))
    expect_equal(ld$Lineage, c("Europe", "NorthAmerica"))
    ## list columns (e.g. segment_color) are not carried over
    expect_false("segment_color" %in% colnames(ld))
    ## viruses without label are filtered out
    expect_false("virus C" %in% ld$label)
})

test_that("prepare_facet_flow_info duplicates flows crossing facets", {
    virus_info <- data.frame(
        id = c("A", "B", "C"),
        Lineage = c("EU", "EU", "NA"),
        stringsAsFactors = FALSE
    )
    flow_info <- build_flow_info(from = c("A", "B"), to = c("B", "C"))

    out <- seqcombo:::prepare_facet_flow_info(virus_info, flow_info, "Lineage")

    expect_equal(nrow(out), 3)
    expect_equal(out$Lineage[out$from == "A"], "EU")
    expect_equal(sort(out$Lineage[out$from == "B"]), c("EU", "NA"))

    ## flows already carrying the facet column are untouched
    flow_info$Lineage <- "EU"
    expect_equal(
        seqcombo:::prepare_facet_flow_info(virus_info, flow_info, "Lineage"),
        flow_info
    )

    ## unknown facet columns are rejected
    expect_error(
        seqcombo:::prepare_facet_flow_info(virus_info, flow_info, "host"),
        "facet column 'host' not available"
    )
})
