# 2-copulas julia

An implementation of bivariate copulas and bivariate distributions in Julia.

This module has the following capabilities:
* Construction of continous bivaraite distributions with Distributions.jl
* Evaluation of cdf and density
* Conditioning (cdfs)
* Sampling
* Bivariate scatter plots, cdf and density plots and contour plots
* Common copula families:
  * π (independence)
  * M (maximal correlation)
  * W (minimal correlation)
  * Gaussian (correlation coefficient r; r = -1 for W, r = 0 for π, r = 1 for M)
  * Frank (s is real; -inf for W, 0 for π, inf for M)
  * Clayton (t>=-1; -1 for W, 0 for π and inf for M)

## Installation

This is not yet a registered julia package. However this package may be installed using the Julia package manager:

```Julia
julia> ]
pkg> add https://github.com/AnderGray/BivariateCopulas.jl.git
```

## Sklar's theorem

A bivariate copula (2-copula) C is a function C:[0,1]<sup>2</sup> -> [0,1] with the following properties: 

<img src="https://imgur.com/8mDL4fr.png" data-canonical-src="https://imgur.com/8mDL4fr.png" width="500" />

Three important 2-copulas are:

<img src="https://imgur.com/bCIWIDX.png" data-canonical-src="https://imgur.com/bCIWIDX.png" width="400" />

with W and M being bounds on all 2-copulas: W ≤ C ≤ M

Copulas are mainly used in dependence modelling, and can be used to construct any continous multivariate distribution function (df) given their univariate marginals. Most probabilistic dependence problems can be reduced to copulas. This is enabled by a theorem from Sklar:


<img src="https://imgur.com/5D1QOif.png" data-canonical-src="https://imgur.com/5D1QOif.png" width="700" />


Moreover, X and Y are stochastically independent if and only if C<sub>XY</sub> = π, maximally correlated iff C<sub>XY</sub> = M, and minimally correlated iff C<sub>XY</sub> = W.

