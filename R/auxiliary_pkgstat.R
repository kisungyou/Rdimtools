#' Show the number of functions for \pkg{Rdimtools}.
#'
#' This function is mainly used for tracking progress for this package.
#'
#' @examples
#' ## run with following command
#' aux.pkgstat()
#'
#' @rdname aux_pkgstat
#' @export
aux.pkgstat <- function(){
  ndo  = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "do."))))
  nest = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "est."))))
  naux = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "aux."))))
  noos = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "oos."))))

  mdo  = paste("*  cat1{do. } manifold learning techniques           : ",ndo,sep="")
  mest = paste("*  cat2{est.} intrinsic dimension estimation methods : ",nest,sep="")
  moos = paste("*  cat3{oos.} out-of-sample projection methods       : ",noos,sep="")
  maux = paste("*  cat4{aux.} auxiliary functions available          : ",naux,sep="")
  print("* Number of functions available in Rdimtools package")
  print(mdo)
  print(mest)
  print(moos)
  print(maux)
}
