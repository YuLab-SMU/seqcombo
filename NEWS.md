# seqcombo 1.35.1

+ improve ggplot2 compatibility by replacing deprecated `aes_()` usage and switching line border sizing to `linewidth`
+ add validation for `flow_info` ids and empty `segment_color` entries, and preserve existing coordinates for nodes not present in `flow_info`
+ add `build_virus_info()` and `build_flow_info()` helpers, `set_layout(..., preserve_x / preserve_y)`, and curved reassortment links
+ add `build_flow_info_from_long()`, segment labels inside virus glyphs, and `build_segment_caption()` for figure annotations
+ add link aesthetics mapped from `flow_info`, a dedicated `layout_timeline()` helper, and faceting support via `facet_by`
+ refresh package metadata, including the package title, description, and `BugReports` URL, and add the `rlang` import
+ migrate the package vignette from `.Rmd` to Quarto `.qmd` and update the vignette build workflow
+ standardize historical `NEWS.md` entries and fix typos for older releases
+ replace the legacy `docs/`-based GitHub Pages flow with deployment from a dedicated `gh-pages` branch via GitHub Actions
+ remove `deploy.sh` and `seqcombo.Rproj`

# seqcombo 1.34.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# seqcombo 1.32.0

+ Bioconductor RELEASE_3_22 (2025-11-01, Sat)

# seqcombo 1.30.0

+ Bioconductor RELEASE_3_21 (2025-04-17, Thu)

# seqcombo 1.28.0

+ Bioconductor RELEASE_3_20 (2024-10-30, Wed)

# seqcombo 1.26.0

+ Bioconductor RELEASE_3_19 (2024-05-15, Wed)

# seqcombo 1.24.0

+ Bioconductor RELEASE_3_18 (2023-10-25, Wed)

# seqcombo 1.18.0

+ Bioconductor RELEASE_3_16 (2022-11-02, Wed)

# seqcombo 1.17.1

+ update docs (2021-12-15, Wed)
+ remove code that was incorporated into ggmsa

# seqcombo 1.16.0

+ Bioconductor RELEASE_3_14

# seqcombo 1.15.1

+ import yulab.utils (2021-08-20, Fri)
+ move `seqdiff` and `simplot` to the ggmsa package

# seqcombo 1.14.0

+ Bioconductor RELEASE_3_13

# seqcombo 1.12.0

+ Bioconductor RELEASE_3_12 (2020-10-28, Wed)

# seqcombo 1.5.1

+ fix R CMD check by importing `dplyr::n` (2019-01-02, Wed)

# seqcombo 1.1.1

+ better simplot implementation (2018-01-09, Fri)

# seqcombo 0.99.11

+ add `geom_genotype` (2017-08-29, Tue)

# seqcombo 0.99.10

+ add `geom_hybrid` (2017-08-17, Thu)

# seqcombo 0.99.9

+ add `hybrid_plot` (2017-06-30, Fri)

# seqcombo 0.0.3

+ add more parameters for plot, `by`, `xlab`, `color`, and `fill`

# seqcombo 0.0.2

+ add vignette

# seqcombo 0.0.1

+ initial version with a plot method for nucleotide differences between two aligned sequences (2016-11-16, Wed)
