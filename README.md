# SPprm
 - An R package for statistical test with power calculation for paired repeated measures (PRM) comparisons
 - Reference
    - Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2018) Statistical test
      with sample size and power calculation for paired repeated
      measures designs of method comparison studies. tech rep.
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
 ## size
 rho = sqrt(s2w+s2b+mu^2)
 PRMtest(X, rho)$p.val
 ## power
 PRMtest(X, rho+0.5)$p.val

```
