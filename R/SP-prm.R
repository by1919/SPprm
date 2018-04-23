## Internal function to compute the log likelihood for paired data Dij under the null value of rho
## @param xpar vector of log(s2w), log(s2b), mu0
## @param Dij  n by m matrix of paired measures
## @param rho  null value of sqrt(s2w+s2b+mu0^2)
llkfun <- function(xpar, Dij, rho=3){
   n = dim(Dij)[1]; m = dim(Dij)[2]
   s2w = exp(xpar[1])
   s2b = exp(xpar[2])
   mu0 = xpar[3]
   E = Dij-mu0
   E1 = rowMeans(E)
   ## Sinv = 1/s2w*diag(m) + (1/(m*s2b+s2w)-1/s2w)/m
   ## det(S): (m*s2b+s2w)*(s2w)^(m-1)
   llk1 = sum(E^2)/s2w + sum(E1^2)*m*(1/(m*s2b+s2w)-1/s2w)
   llk2 = log((m*s2b+s2w)) + (m-1)*xpar[1]
   ##
   llk1+n*llk2
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
## P-value calc with finite sample adjustment
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
#' @return
#' \describe{
#'   \item{p.value}{ test p-values for: RMS test, score Z-test, Wald Z-test }
#' }
#' @export
#' @references
#' Bai,Y., Wu,B. and ... (2018) Statistical test with sample size and power calculation for paired repeated measures designs of method comparison studies.
#' @examples
#' n=20; m=20; s2w=4.44; s2b=1.2; mu0=0.94; rho=3
#' e = matrix(rnorm(n*m), n,m)*sqrt(s2w)
#' u = rnorm(n)*sqrt(s2b)
#' X = mu0 + u + e
#' PRMtest(X,3)
#' PRMtest(X, sqrt(s2w+s2b+mu0^2))
PRMtest <- function(X, rho=3){
  n = dim(X)[1]; m = dim(X)[2]
  ## est
  fn = function(xpar) llkfun(xpar, X, rho)
  hin = function(xpar) Cfun(xpar, rho)
  aaa = nloptr::cobyla(c(log(5),log(5),1), fn, hin=hin)
  xpar = aaa$par
  e2w = exp(xpar[1]);  e2b = exp(xpar[2]);  eu0 = xpar[3]
  ## rhob = sqrt(e2w+e2b+eu0^2)
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
  zpar = newuoa(xpar, fn)$par
  z2w = exp(zpar[1])
  z2b = exp(zpar[2])
  zu0 = zpar[3]
  V = (2*(z2w+m*z2b)^2+2*(m-1)*z2w^2+4*m*(z2w+m*z2b)*zu0^2)/n/m^2
  Z = Zn/sqrt(V)
  pval[3] = 2*pnorm(-abs(Z))
  return(list(p.value=pval, pars=c(s2w=e2w,s2b=e2b,mu=eu0)) )
}

#' Sample size and power calculation for PRM design
#' 
#' @param alpha desired significance level. Default to 0.05
#' @param n  number of subjects
#' @param m  number of repeated measures for each subject
#' @param s2w  within subject variation
#' @param s2b  between subject variation
#' @param rho   null threshold of acceptable RMS value
#' @param rhoa  actual RMS value
#' @return
#' \describe{
#'   \item{xq}{ computed quantile for the RMS test statistic }
#'   \item{size}{ computed actual type I error }
#'   \item{pwr}{ computed power }
#' }
#' @export
#' @references
#' Bai,Y., Wu,B. and ... (2018) Statistical test with sample size and power calculation for paired repeated measures designs of method comparison studies.
PRMsp <- function(alpha=0.05,n=20, m=20,s2w=4.44, s2b=1.82, rho=3, rhoa=2.7){
  mu0 = sqrt(rho^2-s2w-s2b)
  mu1 = sqrt(rhoa^2-s2w-s2b)
  f0 = function(x) Tpval0(x, n,m,s2w,s2b,mu0)/alpha - 1
  xq = uniroot(f0, c(1,1e5))$root
  size = Tpval0(xq, n,m,s2w,s2b,mu0)
  pwr1 = Tpval1(xq, n,m,s2w,s2b,mu1)
  return(list(xq=xq,size=size,pwr=pwr1))
}


Tpval1 <- function(xq, n=20, m=20, s2w=4.44, s2b=1.82, mu0=1.11){
  ## lam = c(s2w+m*s2b, s2w, -xq/(n-1));   h = c(n,n*(m-1), n-1) ## too conservative
  ## lam = c(s2w+m*s2b, s2w, -xq/n/(m-1));   h = c(n,n*(m-1), n*(m-1))
  ## lam = c(s2w+m*s2b, s2w, -xq/(n*m-1));   h = c(n,n*(m-1), n*m-1)
  df1 = n*m-Edf(n,m,s2w,s2b)
  lam = c(s2w+m*s2b, s2w, -xq/df1);   h = c(n,n*(m-1), df1)
  eta0 = n*m*mu0^2/lam[1]
  dta = c(eta0, 0,0)
  1-pchisum(0,lam,h,dta)
}
## Empirical DF for the variance parameters.
Edf <- function(n,m,s2w,s2b){
  X2 = diag(n+1)*m 
  X2[1,] = X2[,1] = m
  X2[1] = n*m
  X2r = X2
  diag(X2r)[-1] = diag(X2r)[-1] + s2w/s2b
  a = solve(X2r,X2)
  return(sum(diag(a)))
}
