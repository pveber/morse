# MORSE: MOdelling tools for Reproduction and Survival data in Ecotoxicology

Companion R package for MOSAIC website
 - need `Rtools 3.3` and last `R` version.

## Versions

### CRAN release version

[![CRAN version](http://www.r-pkg.org/badges/version/morse)](http://cran.rstudio.com/web/packages/morse/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/morse)](http://cran.rstudio.com/web/packages/morse/index.html)

### Build status for development version

[![Build Status](https://api.travis-ci.org/philipperuiz/morse.svg?branch=v2.2.0Complet)](https://api.travis-ci.org/philipperuiz/morse.svg)

## Report a problem

Use the [issue tracker](https://github.com/pveber/morse/issues)
to report a problem.


## Compilation (for command-line users)

- `cd` to source directory
- Build from sources (creates vignette and archive)
  `R CMD build .`
- Check the package
  `R CMD check --as-cran morse_X.X.X.tar.gz`
- Update package description/NAMESPACE
  - under the R interpreter: `roxygen2::roxygenise(".")`
- Generate documentation
  - reference manual: `R CMD Rd2pdf .`
  - vignettes (using the R interpreter):
    `devtools::document(roclets=c('rd', 'collate', 'namespace', 'vignette'))`
- Run unit tests
  - under the R interpreter: `devtools::test()`

## Compilation (with RStudio)

- need `devtools`, `ROxygen2 v5.0.1 or higgher`
- RStudio builder configuration:
  - Project Options :
      enable `Use devtools package...`
      enable `Generate documentation...`
      ROxygen options...`:
        all enable exept `Vignettes` and `Source and binary package build`
  - vignettes folder must contain only 3 files:
      biblio.bib
      modelling.Snw
      tutorial.Rmd
  - update Documentation (no vignette, only NAMESPACE and Rd)
      `Document` or Ctrl + Shift + D
  - check if the two folder inst and build were created.
  - build the source file :
    `More : Build Source Package`

## Install the package from source
R: remove.packages("morse", lib="~/Documents/R/R-X.X.X/library")
R: install.packages("~/Documents/morse_X.X.X.tar.gz", repos = NULL, type = "source")
R: library("morse")
R: ?morse

## Tests
- Build from sources: `Build & Reload`
- Test: `Build > More > Test package` (Ctrl + Shift + T). Note that you must activate the option "Use devtools package functions if available" in `Project Options > Build Tools`.
- Generate the documentation: `Build & Reload or Build > More > Document` (Ctrl + Shift + D)
- Check the package with the stable and the devel version of R:
  - in `Build > More > Configure Build Tools`, add `--as-cran` in `Check Package -- R CMD check additional` options 
  - `Build > Check` (Ctrl + Shift + E)
  

