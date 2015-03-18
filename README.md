r-mosaic
========

Companion R package for MOSAIC website

How to :
========

## Generate documentation : (this action can be run without build or check the package) 
  (with the library roxygen2)

### update or create Rd file and NAMESPACE :
  * R : (in the package directory) `roxygen2::roxygenise(".")`
  * RStudio : Build & Reload or Build > More > Document (Ctrl + Shift + D)

### pdf :
  * cmd (in the parent directory) : `R CMD Rd2pdf mosaic-r`

## Check the package (only when the package was complete because of the tests and vignettes) :

### cmd :
  * 1) build the source : `R CMD build mosaic-r` (it create the vignette and the
       tar.gz)
  * 2) check the package :`R CMD check --as-cran morse_X.X.X.tar.gz`
  
### RStudio :
  * 1) Build & Reload : load the package and update Rd files
  * 2) BUild > Check : build the tar.gz and check (Ctrl + Shift + E)
     (in Build > More > Configure Build Tools :
     add `--as-cran` in Check Package -- R CMD check additional options)

## Test the package :
(with the library testthat)

### cmd :
  * `R CMD check morse_X.X.X.tar.gz`

### RStudio :
  * Build > More > Test package (Ctrl + Shift + T)
