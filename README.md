# BivariateCopulas.jl

An implementation of bivariate copulas and bivariate distributions in Julia.

This module has the following capabilities:
* Construction of continous joint distribution with Distributions.jl
* Evaluations of cdf and density
* Conditioning
* Sampling
* Bivariate scatter plots, cdf and density plots and contour plots
* Common copula families:
  * π (independence)
  * M (maximal correlation)
  * W (minimal correlation)
  * Gaussian (correlation coefficient r; r = -1 for W, r = 0 for π, r = 1 for M)
  * Frank (s is real; -inf for W, 0 for π, inf for M)
  * Clayton (t>-1; -1 for W, 0 for π and inf for M)
