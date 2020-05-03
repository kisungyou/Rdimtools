#' Build a centered kernel matrix K
#'
#' From the celebrated Mercer's Theorem, we know that for a mapping \eqn{\phi}, there exists
#' a kernel function - or, symmetric bilinear form, \eqn{K} such that \deqn{K(x,y) = <\phi(x),\phi(y)>} where \eqn{<,>} is
#' standard inner product. \code{aux.kernelcov} is a collection of 20 such positive definite kernel functions, as
#' well as centering of such kernel since covariance requires a mean to be subtracted and
#' a set of transformed values \eqn{\phi(x_i),i=1,2,\dots,n} are not centered after transformation.
#' Since some kernels require parameters - up to 2, its usage will be listed in arguments section.
#'
#' @details
#' There are 20 kernels supported. Belows are the kernels when given two vectors \eqn{x,y}, \eqn{K(x,y)}
#' \describe{
#' \item{linear}{\eqn{=<x,y>+c}}
#' \item{polynomial}{\eqn{=(<x,y>+c)^d}}
#' \item{gaussian}{\eqn{=exp(-c\|x-y\|^2)}, \eqn{c>0}}
#' \item{laplacian}{\eqn{=exp(-c\|x-y\|)}, \eqn{c>0}}
#' \item{anova}{\eqn{=\sum_k exp(-c(x_k-y_k)^2)^d}, \eqn{c>0,d\ge 1}}
#' \item{sigmoid}{\eqn{=tanh(a<x,y>+b)}}
#' \item{rational quadratic}{\eqn{=1-(\|x-y\|^2)/(\|x-y\|^2+c)}}
#' \item{multiquadric}{\eqn{=\sqrt{\|x-y\|^2 + c^2}}}
#' \item{inverse quadric}{\eqn{=1/(\|x-y\|^2+c^2)}}
#' \item{inverse multiquadric}{\eqn{=1/\sqrt{\|x-y\|^2+c^2}}}
#' \item{circular}{\eqn{=
#' \frac{2}{\pi} arccos(-\frac{\|x-y\|}{c}) - \frac{2}{\pi} \frac{\|x-y\|}{c}\sqrt{1-(\|x-y\|/c)^2}
#' }, \eqn{c>0}}
#' \item{spherical}{\eqn{=
#' 1-1.5\frac{\|x-y\|}{c}+0.5(\|x-y\|/c)^3
#' }, \eqn{c>0}}
#' \item{power/triangular}{\eqn{=-\|x-y\|^d}, \eqn{d\ge 1}}
#' \item{log}{\eqn{=-\log (\|x-y\|^d+1)}}
#' \item{spline}{\eqn{=
#' \prod_i (
#' 1+x_i y_i(1+min(x_i,y_i)) - \frac{x_i + y_i}{2} min(x_i,y_i)^2
#' + \frac{min(x_i,y_i)^3}{3}
#' )
#' }}
#' \item{Cauchy}{\eqn{=\frac{c^2}{c^2+\|x-y\|^2}}}
#' \item{Chi-squared}{\eqn{=\sum_i \frac{2x_i y_i}{x_i+y_i}}}
#' \item{histogram intersection}{\eqn{=\sum_i min(x_i,y_i)}}
#' \item{generalized histogram intersection}{\eqn{=sum_i min(
#' |x_i|^c,|y_i|^d
#' )}}
#' \item{generalized Student-t}{\eqn{=1/(1+\|x-y\|^d)}, \eqn{d\ge 1}}
#' }
#'
#' @param X an \eqn{(n\times p)} data matrix
#' @param ktype a vector containing the type of kernel and parameters involved. Below the usage is
#' consistent with description
#' \describe{
#' \item{linear}{\code{c("linear",c)}}
#' \item{polynomial}{\code{c("polynomial",c,d)}}
#' \item{gaussian}{\code{c("gaussian",c)}}
#' \item{laplacian}{\code{c("laplacian",c)}}
#' \item{anova}{\code{c("anova",c,d)}}
#' \item{sigmoid}{\code{c("sigmoid",a,b)}}
#' \item{rational quadratic}{\code{c("rq",c)}}
#' \item{multiquadric}{\code{c("mq",c)}}
#' \item{inverse quadric}{\code{c("iq",c)}}
#' \item{inverse multiquadric}{\code{c("imq",c)}}
#' \item{circular}{\code{c("circular",c)}}
#' \item{spherical}{\code{c("spherical",c)}}
#' \item{power/triangular}{\code{c("power",d)}}
#' \item{log}{\code{c("log",d)}}
#' \item{spline}{\code{c("spline")}}
#' \item{Cauchy}{\code{c("cauchy",c)}}
#' \item{Chi-squared}{\code{c("chisq")}}
#' \item{histogram intersection}{\code{c("histintx")}}
#' \item{generalized histogram intersection}{\code{c("ghistintx",c,d)}}
#' \item{generalized Student-t}{\code{c("t",d)}}
#' }
#'
#' @return a named list containing \describe{
#' \item{K}{a \eqn{(p\times p)} kernelizd gram matrix.}
#' \item{Kcenter}{a \eqn{(p\times p)} centered version of \code{K}.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate a toy data
#' set.seed(100)
#' X = aux.gensamples(n=100)
#'
#' ## compute a few kernels
#' Klin = aux.kernelcov(X, ktype=c("linear",0))
#' Kgau = aux.kernelcov(X, ktype=c("gaussian",1))
#' Klap = aux.kernelcov(X, ktype=c("laplacian",1))
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(Klin$K, main="kernel=linear")
#' image(Kgau$K, main="kernel=gaussian")
#' image(Klap$K, main="kernel=laplacian")
#' par(opar)
#' }
#'
#'
#' @references Hofmann, T., Scholkopf, B., and Smola, A.J. (2008) \emph{Kernel methods in
#' machine learning}. arXiv:math/0701907.
#' @author Kisung You
#' @rdname aux_kernelcov
#' @export
aux.kernelcov <- function(X,ktype){
  # 6-1. 20 ktype is supported
  if (length(ktype)==1){
    defaultflag = TRUE
  } else {
    defaultflag = FALSE
  }
  if (length(ktype)>3){
    stop("* aux.kernelcov : maximum 2 parameters are used.")
  }
  if (typeof(ktype[1])!="character"){
    stop("* aux.kernelcov : first argument should be a name of kernel chosen.")
  }
  if ((ktype[1]=="circular")&&(ncol(X)!=2)){
    stop("* aux.kernelcov : circular kernel is for 2-d data only.")
  }
  if ((ktype[1]=="spherical")&&(ncol(X)!=3)){
    stop("* aux.kernelcov : spherical kernel is for 3-d data only.")
  }

  kchrname = ktype[1]
  switch(kchrname,
         # 1. "linear" - linear kernel
         linear={
           knumber = as.integer(1)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         #  2. "polynomial" - polynomial kernel
         polynomial={
           knumber = as.integer(2)
           if (defaultflag==TRUE){
             par1 = 1.0
             par2 = 3
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'polynomial' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
           }
         },
         #  3. "gaussian" - gaussian kernel
         gaussian={
           knumber = as.integer(3)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
             if ((par1 <= 0)||(is.na(par1))||(is.infinite(par1))){
               stop("* aux.kernelcov : 'gaussian' scale parameter should be a positive value.")
             }
           }
           par2 = 0.0
         },
         #  4. "laplacian" - laplacian kernel
         laplacian = {
           knumber = as.integer(4)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
             if ((par1 <= 0)||(is.na(par1))||(is.infinite(par1))){
               stop("* aux.kernelcov : 'gaussian' scale parameter should be a positive value.")
             }
           }
           par2 = 0.0
         },
         #  5. "anova" - ANOVA kernel
         anova = {
           knumber = as.integer(5)
           if (defaultflag==TRUE){
             par1 = 1.0
             par2 = 1.0
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'anova' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
             if ((par1<=0)||(par2<1)){
               stop("* aux.kernelcov : 'anova' parameters are invalid.")
             }
           }
         },
         #  6. "sigmoid"
         sigmoid = {
           knumber = as.integer(6)
           if (defaultflag==TRUE){
             par1 = 1/nrow(X)
             par2 = 1.0
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'sigmoid' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
           }
         },
         #  7. "rq" - rational quadratic
         rq = {
           knumber = as.integer(7)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         #  8. "mq" - multiquadric
         mq = {
           knumber = as.integer(8)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         #  9. "iq" - inverse quadric
         iq = {
           knumber = as.integer(9)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         # 10. "imq" - inverse multiquadric
         imq = {
           knumber = as.integer(10)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         # 11. "circular" - circular kernel
         #    use par2 as pi = 3.141592
         circular = {
           knumber = as.integer(11)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           if (par1 <= 0){
             stop("* aux.kernelcov : parameter has to be a positive number.")
           }
           par2 = pi
           if (ncol(X)!=2){
             stop("* aux.kernelcov : 'circular' kernel should be used on 2-dimensional data.")
           }
         },
         #  12. "spherical" - spherical kernel
         #      use par2 as pi = 3.141592
         spherical = {
           knumber = as.integer(12)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           if (par1 <= 0){
             stop("* aux.kernelcov : parameter has to be a positive number.")
           }
           par2 = pi
           if (ncol(X)!=3){
             stop("* aux.kernelcov : 'spherical' kernel should be used on 3-dimensional data.")
           }
         },
         #  13. "power" - power/triangular kernel
         power = {
           knumber = as.integer(13)
           if (defaultflag==TRUE){
             par1 = 2.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1 < 1){
             stop("* aux.kernelcov : 'power' parameter should be >= 1.")
           }
         },
         #  14. "log" - log kernel
         log = {
           knumber = as.integer(14)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1 < 1){
             stop("* aux.kernelcov : 'log' parameter should be >= 1.")
           }
         },
         #  15. "spline" - No Parameter Kernel
         spline = {
           knumber = as.integer(15)
           par1 = 0.0
           par2 = 0.0
         },
         #  16. "cauchy" - cauchy kernel
         cauchy = {
           knumber = as.integer(16)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1 < 0){
             stop("* aux.kernelcov : 'cauchy' parameter should be >= 0.")
           }
         },
         #  17. "chisq" - Chi-Squared kernel
         chisq = {
           knumber = as.integer(17)
           par1 = 0.0
           par2 = 0.0
         },
         #  18. "histintx" - histogram intersection
         histintx = {
           knumber = as.integer(18)
           par1 = 0.0
           par2 = 0.0
           if (any((X<0))){
             stop("* aux.kernelcov : 'histintx' can be used for data with positive values only.")
           }
         },
         #  19. "ghistintx" - generalized histogram intersection
         ghistintx = {
           knumber = as.integer(19)
           if (defaultflag==TRUE){
             par1 = 1.0
             par2 = 1.0
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'ghistintx' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
           }
         },
         #  20. "t" - student-t distribution
         t = {
           knumber = as.integer(20)
           if (defaultflag==TRUE){
             par1 = 2.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1<1){
             stop("* aux.kernelcov : for student t distribution, we need a parameter >= 1.")
           }
         },
         stop("* aux.kernelcov : invalid kernel type. choose one of three or leave it blank.")
  )
  tX = t(X)
  result = aux_kernelcov(tX,knumber,par1,par2)
  return(result);
}
