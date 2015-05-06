# MORSE: MOdelling tools for Reproduction and Survival data in Ecotoxicology

Companion R package for MOSAIC website

### Compilation (for command-line users)

- Build from sources (creates vignette and archive)
  `R CMD build mosaic-r`
- Test (need `Rtools 3.1`, works with `Rtools 3.2` but generate a warning message):
  - generate the documentation
  - in the package directory under the R interpreter: `devtools::test()` 
  - in the package directory under the R interpreter: `devtools::test()`
  
- Generate documentation
  - under the R interpreter: `roxygen2::roxygenise(".")`
  - for a PDF pversion: `R CMD Rd2pdf morse`
- Check the package
  `R CMD check --as-cran morse_X.X.X.tar.gz`

### Compilation (with RStudio)

- Build from sources: `Build & Reload`
- Test: `Build > More > Test package` (Ctrl + Shift + T). Note that you must activate the option "Use devtools package functions if available" in `Project Options > Build Tools`.
- Generate the documentation: `Build & Reload or Build > More > Document` (Ctrl + Shift + D)
- Check the package:
  - in `Build > More > Configure Build Tools`, add `--as-cran` in `Check Package -- R CMD check additional` options 
  - `Build > Check` (Ctrl + Shift + E)

