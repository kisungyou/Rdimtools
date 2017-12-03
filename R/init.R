.pkgenv <- new.env(parent = emptyenv())

.onAttach <- function(libname, pkgname){
  packageStartupMessage("\n * Rdimtools for Dimension Reduction and Estimation maintained by Kisung You (University of Notre Dame).
                         \n * We are continuously developing this package. Please feel free to share your comments, suggestions or bugs to the maintainer.")
}

# .onLoad <- function(libname, pkgname){
#
#   # the below RStudio version control code is adapted from Dirk Eddelbuettel's library(anytime)
#   # and https://github.com/eddelbuettel/anytime/blob/ef8b1e52b80a99e96f46232dfe29180686327887/R/init.R#L49-L52
#   # in light of this helpful discussion:
#   # http://stackoverflow.com/questions/43247649/rcpparmadillo-bigmemory-crashes-windows-rstudio-but-no-other-gui-os-type
#   # mistakes are of course ours
#
#   # bigKRLS requires R 3.3.0 or newer, which is taken care of by the package description
#   # the code below enables an analogous check on RStudio version
#
#   RStudio_outofdate <- if(.Platform$GUI == "RStudio"){
#     threshold <- if(.Platform$OS.type == "unix") "1.0.136" else "1.1.129"
#     !(Sys.getenv("RSTUDIO", unset="0") == "1" &&
#         exists("RStudio.Version") &&
#         eval(parse(text=paste0("RStudio.Version()$version ",  " >= ", "\"", threshold, "\""))))
#   }
#
#   .pkgenv[["RStudio_outofdate"]] <- RStudio_outofdate
# zzzz = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "do."))))
# }


# GITHUB
# git remote add origin git@github.com:kisungyou/Rdimtools.git <- called origin
# git push -u origin master
