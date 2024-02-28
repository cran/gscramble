gscramble
================

<!-- badges: start -->

[![check-standard](https://github.com/eriqande/gscramble/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/eriqande/gscramble/actions/workflows/check-standard.yaml)
[![pkgdown](https://github.com/eriqande/gscramble/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/eriqande/gscramble/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

This is an R package for simulating individuals that are admixed between
different populations according to a pedigree. To do so it uses sampling
without replacement from samples taken from each population.

You can read the documentation for it at:
<https://eriqande.github.io/gscramble/>

You can install in from GitHub:

``` r
remotes::install_github(
  "eriqande/gscramble", 
  build_opts = c("--no-resave-data"), 
  build_vignettes = TRUE, 
  build_manual = TRUE,
  quiet = FALSE
)
```

Once that is installed you can read the vignettes:

``` r
vignette("gscramble-tutorial")  # information on input data objects
vignette("about-createGSP")  # a function to make simple genomic simulation pedigrees
vignette("permutation-options")  # about the different options for permuting
```
