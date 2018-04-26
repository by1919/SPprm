# SPprm
 - An R package for statistical test with power calculation for paired repeated measures (PRM) comparisons
 - Reference
    - Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2018) Statistical test
      with sample size and power calculation for paired repeated
      measures designs of method comparison studies. *tech rep*.
 - Sample R codes
```R
 ## install the package
 devtools::install_github('baolinwu/SPprm')
 library(SPprm)
 ## simulate PRM data
 n = 20;  m = 10
 s2w = s2b = 4;  mu = 1
 e = matrix(rnorm(n*m), n,m)*sqrt(s2w)
 u = rnorm(n)*sqrt(s2b)
 X = mu + u + e
 ## QMS test
 rho = sqrt(s2w+s2b+mu^2)
 PRMtest(X, rho)$p.val
 PRMtest(X, rho+0.75)$p.val
 ## power calculation
 ## analytical calc
 PRMap(alpha=0.05,n, m, s2w, s2b, mu, rho+0.75)$pwr
 ## Monte Carlo sim
 PRMmcp(alpha=0.05,n, m, s2w, s2b, mu, rho+0.75)$pwr
 ## sample size calc
 aa = PRMas(pwr=0.8,alpha=0.05, m,s2w, s2b, mu, rho+0.75)
 aa
 ## PRMap(alpha=0.05,n=aa$n, m, s2w, s2b, mu, rho+0.75)$pwr
 PRMmcp(alpha=0.05,n=aa$n, m, s2w, s2b, mu, rho+0.75)$pwr

```
