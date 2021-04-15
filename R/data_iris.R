#' 'Iris' data
#'
#' This is the identical dataset as original \code{iris} data where numeric values of
#' \code{Sepal.Length}, \code{Sepal.Width}, \code{Petal.Length}, \code{Petal.Width}
#' measured in centimeters are given for 50 flowers from each of 3 species of iris.
#'
#' @usage data(iris)
#'
#' @examples
#' \donttest{
#' # load the data
#' data(iris)
#'
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' plot(iris[,1:4])
#' par(opar)
#' }
#'
#' @format a data.frame containing\describe{
#' \item{Sepal.Length}{sepal length}
#' \item{Sepal.Width}{sepal width}
#' \item{Petal.Length}{petal length}
#' \item{Petal.Width}{petal width}
#' \item{Species}{(factor) one of 'setosa','versicolor', and 'virginica'.}
#' }
#'
#' @concept data
"iris"
