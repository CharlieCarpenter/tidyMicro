
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyMicro

<!-- badges: start -->

<!-- badges: end -->

The goal of tidyMicro is to provide a cohesive pipeline for microbiome
data analysis and visualization.

## Installation

You can install the latest version of tidyMicro from this GitHub using

``` r
devtools::install_github("CharlieCarpenter/tidyMicro", 
                         build_vignettes = TRUE) ## Recommended
```

Or you can install the latest release of tidyMicro from
[CRAN](https://CRAN.R-project.org) using

``` r
install.packages("tidyMicro")
```

## Main Features

  - Merge any number OTU tables with clinical data into one tidy data
    set

  - Data visualization for exploration

  - Calculate, analyze, and visualize diversity measures

  - Fit negative binomial, beta binomial, rank sum, and presence absence
    models

  - Create stacked bar charts based on negative binomial and beta
    binomial model estimates

  - Longitudinal PCA plots
