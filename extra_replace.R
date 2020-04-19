# after pkgdown, let's set all the families

rm(list=ls())
library(xfun)
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# step 1. list all files in R directory
Rfiles = list.files("R")
Rpaths = paste0(getwd(),"/",Rfiles)

# step 2. do "linear_" filenames
Rlinear = paste0(getwd(),"/R/",list.files("R",pattern="^linear_"))
for (i in 1:length(Rlinear)){
  val = grep("rdname linear_",readLines(Rlinear[i]))
  if (length(val)>0){
    gsub_file(file=Rlinear[i],"@export", "@family linear_methods \n#' @export")
  }
}

# step 3. do "nonlinear_" functions
Rfiles = paste0(getwd(),"/R/",list.files("R",pattern="^nonlinear_"))
for (i in 1:length(Rfiles)){
  val = grep("rdname nonlinear_",readLines(Rfiles[i]))
  if (length(val)>0){
    gsub_file(file=Rfiles[i],"@export", "@family nonlinear_methods \n#' @export")
  }
}


# Correction --------------------------------------------------------------
# If @family is used, all the functions are cross-referenced.
# In order to avoid this, use @concept rather for better grouping.
# It's not automatically added but still way better.

Rfiles = paste0(getwd(),"/R/",list.files("R"))
for (i in 1:length(Rfiles)){
  val = grep("@family",readLines(Rfiles[i]))
  if (length(val)>0){
    xfun::gsub_file(file=Rfiles[i],"@family", "@concept")
  }
}


# Correction 2. figure ----------------------------------------------------
# mfrow=c(3,1) -> mfrow=c(1,3)

Rfiles = paste0(getwd(),"/R/",list.files("R"))
for (i in 1:length(Rfiles)){
  val = grep("mfrow=c(3,1)",readLines(Rfiles[i]))
  if (length(val)>0){
    xfun::gsub_file(file=Rfiles[i],"mfrow=c(3,1)", "mfrow=c(1,3)")
  }
}
