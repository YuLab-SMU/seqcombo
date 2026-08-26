test_that("seqcombo_host_palette returns the full built-in set by default", {
    pal <- seqcombo_host_palette()

    expect_named(pal, c("Human", "Avian", "Swine", "Equine", "Bat"))
    expect_true(all(startsWith(pal, "#")))
})

test_that("seqcombo_host_palette matches labels case-insensitively", {
    pal <- seqcombo_host_palette(c("avian", "Human"))

    expect_identical(
        unname(pal),
        unname(seqcombo_host_palette()[c("Avian", "Human")])
    )
    expect_identical(names(pal), c("avian", "Human"))
})

test_that("seqcombo_host_palette assigns deterministic colors to unknown levels", {
    levels <- c("Human", "platypus")
    first <- seqcombo_host_palette(levels)
    second <- seqcombo_host_palette(rev(levels))

    expect_false(is.na(first[["platypus"]]))
    expect_identical(first[["platypus"]], second[["platypus"]])
    expect_identical(first[["Human"]], seqcombo_host_palette()[["Human"]])
})

test_that("apply_seqcombo_palette maps every label to a color", {
    hosts <- c("Avian", NA, "swine", "Swine")

    cols <- apply_seqcombo_palette(hosts)

    expect_null(names(cols))
    expect_length(cols, length(hosts))
    expect_identical(cols[3], cols[4])
    expect_identical(cols[1], seqcombo_host_palette()[["Avian"]])
    expect_identical(cols[2], "#BEBEBE")
})

test_that("scale_seqcombo_host supplies requested manual scales", {
    scales_list <- scale_seqcombo_host()

    expect_length(scales_list, 2)
    expect_s3_class(scales_list[[1]], "ScaleDiscrete")
    expect_s3_class(scales_list[[2]], "ScaleDiscrete")

    colour_only <- scale_seqcombo_host(aesthetics = "color")
    expect_length(colour_only, 1)
})

test_that("scale_seqcombo_host attaches to an autoplot output", {
    data <- example_seqcombo_data()
    p <- autoplot(data, v_color = ~Host, v_fill = ~Host) +
        scale_seqcombo_host()

    expect_s3_class(p, "ggplot")
})
