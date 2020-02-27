#' Intrinsic Dimensionality Estimation with DANCo
#'
#'
#'
#' @examples
#' \dontrun{
#' ## create 3 datasets of intrinsic dimension 2.
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="saddle")
#'
#' ## acquire an estimate for intrinsic dimension
#' out1 = est.danco(X1, k=10)
#' out2 = est.danco(X2, k=10)
#' out3 = est.danco(X3, k=10)
#'
#' ## print the results
#' sprintf("* est.danco : estimated dimension for 'swiss'  data is %.2f.",out1$estdim)
#' sprintf("* est.danco : estimated dimension for 'ribbon' data is %.2f.",out2$estdim)
#' sprintf("* est.danco : estimated dimension for 'saddle' data is %.2f.",out3$estdim)
#' }
#'
#' @references
#' \insertRef{ceruti_danco_2014}{Rdimtools}
#'
#' @export
est.danco <- function(X, k=5){
  ##########################################################################
  ## preprocessing
  aux.typecheck(X)
  N = nrow(X)
  D = ncol(X)
  k = round(k)

  ##########################################################################
  ## computation part 1
  #   index of (k+1) nearest neighbors
  dX = as.matrix(stats::dist(X))
  vec.topK1 <- list()
  for (n in 1:N){
    tgt = as.vector(dX[n,])
    vec.topK1[[n]] = order(tgt)[2:(k+2)]
  }
  #   Lombardi's procedure
  data.dML = danco_part1(dX, vec.topK1, N, D, k)

  ##########################################################################
  ## computation part 2
  muvec = danco_part2(X, vec.topK1, N, D, k)
  data.vvv = as.double(muvec[1])
  data.tau = as.double(muvec[2])

  ##########################################################################
  ## computation part 3 : bootstrapping
  boot.dML = rep(0,D)
  boot.vvv = rep(0,D)
  boot.tau = rep(0,D)
  for (d in 2:D){
    Y  = matrix(stats::rnorm(N*d),ncol=d)
    Y  = (Y/base::sqrt(base::rowSums(Y^2)))*(stats::runif(1, min=1e-8, max=1-(1e-8))^(1/d))
    dY = as.matrix(stats::dist(Y))
    vec.topK1 <- list()
    for (n in 1:N){
      tgt = as.vector(dY[n,])
      vec.topK1[[n]] = order(tgt)[2:(k+2)]
    }
    boot.dML[d] = danco_part1(dY, vec.topK1, N, D, k)

    muvec = danco_part2(Y, vec.topK1, N, D, k)
    boot.vvv[d] = as.double(muvec[1])
    boot.tau[d] = as.double(muvec[2])
  }

  ##########################################################################
  ## computation 4 : compute all distances
  d.cost = rep(0,D)
  for (d in 2:D){
    d.cost[d] = danco_cost(k, data.dML, boot.dML[d], data.vvv, data.tau, boot.vvv[d], boot.tau[d])
  }
  d.fin = base::which.max(d.cost)
  if (length(d.fin) > 1){
    d.fin = base::sample(d.fin, 1)
  }

  ##########################################################################
  ## Return the results
  result = list()
  result$estdim = d.fin
  return(result)
}


# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
danco_part1 <- function(dmat, vec.topK1, N, D, k){
  vec.rho = rep(0,N)
  for (n in 1:N){
    tgt = dmat[n,vec.topK1[[n]]]
    vec.rho[n] = min(tgt/(max(tgt)))
  }
  ## Algorithm 1 : Simply Searching Over
  # val.lomb = rep(0,D)
  # for (d in 1:D){
  #   term1 = N*log(k*d)
  #   term2 = (d-1)*sum(log(vec.rho))
  #   term3 = (k-1)*sum(log(1-(vec.rho^d)))
  #   val.lomb[d] = term1+term2+term3
  # }
  # idmax = which.max(val.lomb)
  # if (length(idmax) > 1){
  #   return(base::sample(idmax, 1))
  # } else {
  #   return(idmax)
  # }
  ## Algorithm 2 : Optimization
  funmax <- function(d){
    term1 = N*log(k*d)
    term2 = (d-1)*sum(log(vec.rho))
    term3 = (k-1)*sum(log(1-(vec.rho^d)))
    return(term1+term2+term3)
  }
  return(as.double(stats::optimise(funmax, lower=1, upper=D, maximum=TRUE)$maximum))
}
#' @keywords internal
#' @noRd
danco_part2 <- function(X, vec.topK1, N, D, k){
  thetas = c() # let's stack them as row vectors
  for (n in 1:N){
    tgtidx   = vec.topK1[[n]][1:k]
    thetavec = as.vector(danco_part2_theta(X[tgtidx,]))
    thetas   = base::rbind(thetas, thetavec)
  }

  vec.vvv = rep(0,N)
  vec.eta = rep(0,N)
  vec.tau = rep(0,N)
  for (n in 1:N){
    tgt = as.vector(thetas[n,])
    tgt.sin = base::sin(tgt)
    tgt.cos = base::cos(tgt)
    eta     = sqrt((sum(tgt.cos)/N)^2 + (sum(tgt.sin)/N)^2)
    vec.vvv[n] = base::atan2(base::sum(tgt.sin), base::sum(tgt.cos))
    vec.eta[n] = eta
    if (eta < 0.53){
      vec.tau[n] = 2*eta + (eta^3) + (5*(eta^5))/6
    } else if (eta >= 0.85){
      vec.tau[n] = 1/((eta^3) - 4*(eta^2) + (3*eta))
    } else {
      vec.tau[n] = -0.4 + 1.39*eta + 0.43/(1-eta)
    }
  }

  mu.vvv = base::atan2(base::sum(base::sin(vec.vvv)), base::sum(base::cos(vec.vvv)))
  mu.tau = (sum(vec.tau)/N)
  return(c(mu.vvv, mu.tau))
}
#' @keywords internal
#' @noRd
danco_part2_theta <- function(X){
  k = base::nrow(X)
  norms  = rep(0,k)
  for (i in 1:k){
    tgt = as.vector(X[i,])
    norms[i] = base::sqrt(base::sum(tgt^2))
  }
  output = array(0,c(k,k))
  for (i in 1:(k-1)){
    tgti = as.vector(X[i,])
    nrmi = norms[i]
    for (j in (i+1):k){
      tgtj = as.vector(X[j,])
      nrmj = norms[j]
      output[i,j] <- output[j,i] <- base::acos(base::sum(tgti*tgtj)/(nrmi*nrmj))
    }
  }
  return(as.vector(output[base::upper.tri(output)]))
}
#' @keywords internal
#' @noRd
danco_cost <- function(k, dML, ddML, muvvv, mutau, dvvv, dtau){
  # harmonic summation
  vecHk = 1/(1:k)
  Hk  = base::sum(vecHk)
  Hk1 = base::sum(vecHk[1:(k-1)])

  # sub 1
  term1 = Hk*(ddML/dML) - 1 - Hk1 - log(ddML) + log(dML)
  term2 = 0
  factk = base::factorial(k)
  for (i in 0:k){
    term2 = term2 + ((-1)^i)*(factk/(factorial(i)*factorial(k-i)))*base::digamma(1 + i*(dML/ddML))
  }
  sub1 = term1 - (k-1)*term2

  # sub 2
  v1 = muvvv
  v2 = dvvv
  tau1 = mutau
  tau2 = dtau

  term1 = log(besselI(tau2, 0)) - log(besselI(tau1, 0))
  term2 = (intbessel1(tau1) - intbessel1(-tau1))/(2*besselI(tau1, 0))
  term3 = tau1 - tau2*cos(v2-v1)
  sub2  = term1 + (term2*term3)

  # return
  return(sub1+sub2)
}
#' @keywords internal
#' @noRd
intbessel1 <- function(v){
  myfun <- function(theta){
    return(exp(v*cos(theta))*cos(theta))
  }
  return(as.double(stats::integrate(myfun, lower=0, upper=pi)$value)/pi)
}
