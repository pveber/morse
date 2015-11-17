# MORSE: MOdelling tools for Reproduction and Survival data in Ecotoxicology

Companion R package for MOSAIC website
 - need `Rtools 3.3` and last `R` version.

### Compilation (for command-line users (Don't build the vignette
### correctly for now ==> use RStudio only)) 

- Build from sources (creates vignette and archive)
  `R CMD build mosaic-r`
- need `Rtools 3.3`:
  - generate the documentation
  - in the package directory under the R interpreter: `devtools::test()`
  
- Generate documentation
  - under the R interpreter: `roxygen2::roxygenise(".")`
  - for a PDF pversion: `R CMD Rd2pdf morse`
- Check the package
  `R CMD check --as-cran morse_X.X.X.tar.gz`

### Compilation (with RStudio)

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

### Install the package from source
R: remove.packages("morse", lib="~/Documents/R/R-X.X.X/library")
R: install.packages("~/Documents/morse_X.X.X.tar.gz", repos = NULL, type = "source")
R: library("morse")
R: ?morse

### Tests
- Build from sources: `Build & Reload`
- Test: `Build > More > Test package` (Ctrl + Shift + T). Note that you must activate the option "Use devtools package functions if available" in `Project Options > Build Tools`.
- Generate the documentation: `Build & Reload or Build > More > Document` (Ctrl + Shift + D)
- Check the package with the stable and the devel version of R:
  - in `Build > More > Configure Build Tools`, add `--as-cran` in `Check Package -- R CMD check additional` options 
  - `Build > Check` (Ctrl + Shift + E)
  

