## Internal function to compute the log likelihood for paired data Dij under the null value of rho
## @param xpar vector of log(s2w), log(s2b), mu0
## @param Dij  n by m matrix of paired measures
## @param rho  null value of sqrt(s2w+s2b+mu0^2)
llkfun <- function(xpar, Dij, REML=TRUE){
   n = dim(Dij)[1]; m = dim(Dij)[2]
   s2w = exp(xpar[1])
   s2b = exp(xpar[2])
   mu0 = xpar[3]
   ## 
   D1 = rowMeans(Dij)
   E = Dij-D1
   s2s = m*s2b+s2w
   llk0 = log(s2s) + (m-1)*xpar[1]
   llk1 = sum(E^2)/s2w + m*sum((D1-mu0)^2)/s2s
   allk = 0
   if(REML) allk = -log(s2s)
   n*llk0+llk1+allk
}
## Internal constraint function
## @param xpar vector of log(s2w), log(s2b), mu0
## @param rho  null value of sqrt(s2w+s2b+mu0^2)
Cfun <- function(xpar, rho=3){
   s2w = exp(xpar[1])
   s2b = exp(xpar[2])
   mu0 = xpar[3]
   s2w+s2b+mu0^2-rho^2
}
## P-value calc
Tpval0 <- function(xq, n=20, m=20, s2w=4.44, s2b=1.82, mu0=1.11){
  lam = c(s2w+m*s2b, s2w);   h = c(n,n*(m-1))
  eta0 = n*m*mu0^2/lam[1]
  dta = c(eta0, 0)
  1-pchisum(xq,lam,h,dta)
} 


#' Conduct estimation and RMS test for paired repeated measures (PRM) designs of study
#'
#' The problem can be reduced to a mixed effects model estimation under the RMS constraint. Due to the special structure of the
#' mixed effects model for the PRM designs, we can compute the likelihood analytically and hence efficiently calculate
#' MLE. See the reference of Bai et. al (2018). The score and Wald Z-tests are also implemented. 
#'
#' @param X n by m matrix of paired measures
#' @param rho  null threshold of acceptable RMS value
#' @param REML using REML instead of MLE. Default to TRUE.
#' @return
#' \describe{
#'   \item{p.value}{ test p-values for: RMS test, score Z-test, Wald Z-test }
#'   \item{pars}{ estimated null parameter values }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2018) Statistical test with sample size and power calculation for paired repeated measures designs of method comparison studies.
#' @examples
#' n=20; m=20; s2w=4.44; s2b=1.2; mu0=0.94
#' e = matrix(rnorm(n*m), n,m)*sqrt(s2w)
#' u = rnorm(n)*sqrt(s2b)
#' X = mu0 + u + e
#' PRMtest(X,3)
#' PRMtest(X, sqrt(s2w+s2b+mu0^2))
PRMtest <- function(X, rho=3, REML=TRUE){
  n = dim(X)[1]; m = dim(X)[2]
  ## est
  fn = function(xpar) llkfun(xpar, X, REML)
  hin = function(xpar) Cfun(xpar, rho)
  aaa = nloptr::cobyla(c(log(5),log(5),1), fn, hin=hin)
  xpar = aaa$par
  e2w = exp(xpar[1]);  e2b = exp(xpar[2]);  eu0 = xpar[3]
  ## test
  pval = rep(NA,3)
  xq = sum(X^2)
  pval[1] = Tpval0(xq, n,m,e2w,e2b,eu0)
  ## score Z
  Zn = xq/n/m - rho^2
  V0 = (2*(e2w+m*e2b)^2+2*(m-1)*e2w^2+4*m*(e2w+m*e2b)*eu0^2)/n/m^2
  Z = Zn/sqrt(V0)
  pval[2] = 2*pnorm(-abs(Z))
  ## Wald Z
  zpar = nloptr::newuoa(xpar, fn)$par
  z2w = exp(zpar[1])
  z2b = exp(zpar[2])
  zu0 = zpar[3]
  V = (2*(z2w+m*z2b)^2+2*(m-1)*z2w^2+4*m*(z2w+m*z2b)*zu0^2)/n/m^2
  Z = Zn/sqrt(V)
  pval[3] = 2*pnorm(-abs(Z))
  names(pval) = c('RMS', 'Z-score', 'Z-Wald')
  return(list(p.value=pval, pars=c(s2w=e2w,s2b=e2b,mu=eu0)) )
}

#' Analytical power calculation for PRM design problem
#'
#' Due to its nature of composite null hypothesis, it is not trivial to derive the analytical power calculation. We develop an average likelihood
#' based approach.
#' 
#' @param alpha desired significance level. Default to 0.05
#' @param n  number of subjects
#' @param m  number of repeated measures for each subject
#' @param s2w  within subject variation
#' @param s2b  between subject variation
#' @param mu   mean measure difference
#' @param rho   null threshold of acceptable RMS value
#' @param REML using REML instead of MLE. Default to TRUE.
#' @return
#' \describe{
#'   \item{xq}{ computed quantile for the RMS test statistic }
#'   \item{size}{ computed actual type I error }
#'   \item{pwr}{ computed power }
#'   \item{par0}{ estimated null parameter values }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2018) Statistical test with sample size and power calculation for paired repeated measures designs of method comparison studies.
PRMap <- function(alpha=0.05,n=20, m=20,s2w=4.44, s2b=1.82, mu=1.01, rho=3, REML=TRUE){
  ## est
  fn = function(xpar){
    e2w = exp(xpar[1])
    e2b = exp(xpar[2])
    eu0 = xpar[3]
    e2s = m*e2b+e2w
    llk0 = log(e2s) + (m-1)*xpar[1]
    llk1 = (m-1)*s2w/e2w + (m*(eu0-mu)^2 + m*s2b+s2w)/e2s
    allk = 0
    if(REML) allk = -log(e2s)/n
    llk0+llk1+allk
  }
  ## est
  hin = function(xpar) Cfun(xpar, rho)
  aaa = nloptr::cobyla(c(log(5),log(5),1), fn, hin=hin)
  xpar = aaa$par
  e2w = exp(xpar[1]);  e2b = exp(xpar[2]);  eu0 = xpar[3]
  ##  
  A0 = round(1/alpha)
  f0 = function(x) ( Tpval0(x, n,m,e2w,e2b,eu0)/alpha - 1 )*A0
  xq = uniroot(f0, c(1,1e5))$root
  size = Tpval0(xq, n,m,e2w,e2b,eu0)
  pwr1 = Tpval0(xq, n,m,s2w,s2b,mu)
  return( list(xq=xq,size=size,pwr=pwr1, par0=c(s2w=e2w,s2b=e2b,mu=eu0)) )
}


#' Monte Carlo power calculation for PRM design problem
#'
#' Due to its nature of composite null hypothesis, it is not trivial to derive the analytical power calculation. 
#' But Monte Carlo simulation offers a straightforward approach.
#' 
#' @param alpha desired significance level. Default to 0.05
#' @param n  number of subjects
#' @param m  number of repeated measures for each subject
#' @param s2w  within subject variation
#' @param s2b  between subject variation
#' @param mu   mean measure difference
#' @param rho   null threshold of acceptable RMS value
#' @param REML using REML instead of MLE. Default to TRUE.
#' @param B   total number of Monte Carlo simulations. Default to 1000.
#' @return
#' \describe{
#'   \item{pwr}{ computed power }
#'   \item{p.value}{ p-values from Monte Carlo simulations }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2018) Statistical test with sample size and power calculation for paired repeated measures designs of method comparison studies.
PRMmcp <- function(alpha=0.05,n=20, m=20,s2w=4.44, s2b=1.82, mu=1.01, rho=3, REML=TRUE, B=1e3){
  pval = rep(1, B)
  for(it in 1:B){
    e = matrix(rnorm(n*m), n,m)*sqrt(s2w)
    u = rnorm(n)*sqrt(s2b)
    X = mu + u + e
    pval[it] = PRMtest(X, rho, REML)$p.value[1]
  }
  pwr = mean(pval<=alpha)
  return( list( pwr=pwr, p.value=pval) )
}

