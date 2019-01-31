#' Conduct estimation and RMS test for one-way random-effects ANOVA model
#'
#' We are conducting hypothesis test for a composite parameter, the RMS, defined as \eqn{\sqrt{\mu^2+\sigma_b^2+\sigma_w^2}:=\sqrt{\rho} }, where
#' \eqn{\mu} is the overall mean, and \eqn{(\sigma_b^2,\sigma_w^2)} are the between/within-subject variances in the
#'  one-way random-effects ANOVA model, \eqn{y_{ij}=\mu+u_i+\epsilon_{ij}}, where \eqn{u_i\sim N(0,\sigma_b^2) } and
#'  \eqn{\epsilon_{ij}\sim N(0,\sigma_w^2) }. We want to test \eqn{H_0: \rho\ge \rho_0}.
#' We implement a parametric Bootstrap based test with ``exact'' p-value calculation, voiding the need for Bootstrap Monte Carlo simulation.
#' See the reference of Bai et. al (2018). The score and Wald Z-tests, both large-sample normal approximation tests, are also implemented. 
#'
#' @param Y vector of outcomes
#' @param subj subject id (factors). Observations with the sam id are coming from the same individual.
#' @param rho  null threshold of acceptable squared RMS value.
#' @param REML using REML instead of MLE. Default to TRUE.
#' @return
#' \describe{
#'   \item{p.value}{ test p-values for: QMS test, score Z-test, Wald Z-test }
#'   \item{pars0}{ estimated null parameter values }
#'   \item{pars}{ estimated MLE parameter values }
#' }
#' @export
#' @references
#' Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2018) Statistical test with sample size and power calculation for paired repeated measures designs of method comparison studies.
#' @examples
#' s2w=1.4^2; s2b=1.7^2; mu0=-0.4
#' ng = c(10,2,10,10,5,7,9,10)
#' A = rep(1:8, times=ng)
#' Y = mu0 + (rnorm(8)*sqrt(s2b))[A] + rnorm(sum(ng))*sqrt(s2w)
#' RMSt(Y,A)
RMSt <- function(Y, subj, rho=9, REML=TRUE){
  ng = as.vector( table(subj) ); N = sum(ng); K = length(ng)
  mus = tapply(Y, subj, mean)
  sse = sum( tapply(Y, subj, function(yi) sum((yi-mean(yi))^2) ) )
  ## est
  lfn = function(xpar){
    mu = xpar[1]; s2b = exp(xpar[2]); s2w = exp(xpar[3])
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*xpar[3] + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu)^2 )
    ans = ll1 + REML*log( sum(ng/(s2w+ng*s2b)) )
  }
  cfn = function(xpar){
    mu = xpar[1]; s2b = exp(xpar[2]); s2w = exp(xpar[3])
    mu^2+s2b+s2w - rho
  }
  xpar = nloptr::cobyla(c(mean(Y),log(rho/2),log(rho/2)), lfn, hin=cfn)$par
  mu = xpar[1];  s2b = exp(xpar[2]);  s2w = exp(xpar[3])
  ## RMS test
  lam = h = dta = NULL
  for(i in 1:K){
    lam = c(lam, s2w+ng[i]*s2b, s2w);   h = c(h, 1,ng[i]-1)
    eta0 = ng[i]*mu^2/(s2w+ng[i]*s2b)
    dta = c(dta,eta0, 0)
  }
  pvalt = 1-pchisum(sum(Y^2),lam,h,dta)
  ## score Z
  tau2 = ( 2*(s2w+ng*s2b)^2+2*(ng-1)*s2w^2 + 4*ng*(s2w+ng*s2b)*mu^2 )/ng^2
  Zs = (mean(Y^2)-rho)/sqrt(sum(tau2))
  pvals = pnorm(Zs)
  ## Wald Z
  xpar = nloptr::newuoa(c(0,0,0), lfn)$par
  mu1 = xpar[1]; s2b1 = exp(xpar[2]); s2w1 = exp(xpar[3])
  tau2 = ( 2*(s2w1+ng*s2b1)^2+2*(ng-1)*s2w1^2 + 4*ng*(s2w1+ng*s2b1)*mu1^2 )/ng^2
  Zw = (mean(Y^2)-rho)/sqrt(sum(tau2))
  pvalw = pnorm(Zw)
  ##
  pval = c(pvalt,pvals,pvalw)
  names(pval) = c('QMS', 'Z-score', 'Z-Wald')
  return( list(p.value=pval, pars0=c(s2w=s2w,s2b=s2b,mu=mu), pars=c(s2w=s2w1,s2b=s2b1,mu=mu1)) )
}

