seqcombo: sequence recombination visualization
===========================================================================================================================

## Installation

```r
repo = c("https://cran.rstudio.com",
         "https://bioconductor.org/packages/release/bioc",
         "https://guangchuangyu.github.io/drat")
install.packages("seqcombo", repo=repo)
```

## Example

```r
library(seqcombo)
fas <- list.files(system.file("examples","GVariation", package="seqcombo"),
                  pattern="fas", full.names=TRUE)
x <- lapply(fas, seqdiff)
plts <- lapply(x, plot)
plot_grid(plotlist=plts, ncol=1, labels=LETTERS[1:3])
```

![](https://raw.githubusercontent.com/GuangchuangYu/seqcombo/master/inst/figures/GVariation.png)

## Reference

+ <http://link.springer.com/article/10.1007/s11540-015-9307-3>
