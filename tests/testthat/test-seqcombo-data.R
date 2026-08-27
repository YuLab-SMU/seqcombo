test_that("as_seqcombo_data and check_seqcombo_data validate plotting inputs", {
    data <- example_seqcombo_data("basic")

    object <- as_seqcombo_data(data$virus_info, data$flow_info)

    expect_s3_class(object, "seqcombo_data")
    expect_true(check_seqcombo_data(object))
})

test_that("check_seqcombo_data reports missing coordinates when required", {
    virus_info <- data.frame(
        id = c("avian", "human"),
        segment_color = I(list(
            c("purple", "purple"),
            c("red", "red")
        )),
        stringsAsFactors = FALSE
    )

    expect_error(
        check_seqcombo_data(virus_info, require_coordinates = TRUE),
        "'x' and 'y' columns are required in 'virus_info'"
    )
})

test_that("autoplot.seqcombo_data returns a ggplot for flow and genotype workflows", {
    data <- example_seqcombo_data("basic")
    flow_plot <- autoplot(data)

    expect_s3_class(flow_plot, "ggplot")

    genotype_only <- as_seqcombo_data(data$virus_info, NULL)
    genotype_plot <- autoplot(genotype_only)

    expect_s3_class(genotype_plot, "ggplot")
})

test_that("example_seqcombo_data provides a long-format teaching dataset", {
    data <- example_seqcombo_data("long")

    expect_s3_class(data, "seqcombo_data")
    expect_true(check_seqcombo_data(data))

    tables <- attr(data, "long_tables")
    expect_s3_class(tables$segments, "data.frame")
    expect_s3_class(tables$flows, "data.frame")
    expect_setequal(
        unique(tables$segments$id),
        unique(data$virus_info$id)
    )
})

test_that("example_seqcombo_data provides a genotype-only dataset", {
    data <- example_seqcombo_data("genotype")

    expect_null(data$flow_info)
    expect_true(check_seqcombo_data(data))
})
