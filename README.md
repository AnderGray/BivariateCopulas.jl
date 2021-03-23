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


Where F<sub>X</sub> and F<sub>Y</sub> are the cumulative distribution functions (cdf) of X and Y. Moreover, X and Y are stochastically independent if and only if C<sub>XY</sub> = π, maximally correlated iff C<sub>XY</sub> = M, and minimally correlated iff C<sub>XY</sub> = W.

## Usage

### Constructing copulas

```Julia
julia> a = Pi()
Copula ~ π()

julia> b = M()
Copula ~ M()

julia> c = W()
Copula ~ W()

julia> d = Gaussian(0.8)
Copula ~ Gau(r=0.8)

julia> e = Frank(10)
Copula ~ F(s=10.0)

julia> g = Clayton(-0.2)
Copula ~ Cla(t=-0.2)
```
### CDF and Density

```Julia
julia> cop = Gaussian(0.8);

julia> cop(0.2, 0.4) # cdf eval
0-dimensional Array{Float64,0}:
0.17890296062944128

julia> cdf(cop, 0.2, 0.4) # same as above 
0-dimensional Array{Float64,0}:
0.17890296062944128

julia> cop(0.2:0.1:0.8, 0.2:0.1:0.8) # ask for matrix
7×7 Array{Float64,2}:
 0.129034  0.16003   0.178903  0.189859  0.195758  0.198579  0.199683
 0.16003   0.21118   0.247245  0.271358  0.286343  0.294693  0.298579
 0.178903  0.247245  0.300933  0.340944  0.368789  0.386343  0.395758
 0.189859  0.271358  0.340944  0.397584  0.440944  0.471358  0.489859
 0.195758  0.286343  0.368789  0.440944  0.500933  0.547245  0.578903
 0.198579  0.294693  0.386343  0.471358  0.547245  0.61118   0.66003
 0.199683  0.298579  0.395758  0.489859  0.578903  0.66003   0.729034

julia> density(cop, 0.2, 0.4) # density eval at a point
1.347233996821345

julia> density(cop, 0.2:0.1:0.8, 0.2:0.1:0.8) # ask for matrix
7×7 Array{Float64,2}:
 2.28311    1.8543    1.34723   0.888161  0.522393  0.26094   0.0981165
 1.8543     1.88328   1.65627   1.30532   0.917827  0.554958  0.26094
 1.34723    1.65627   1.71488   1.57426   1.28937   0.917827  0.522393
 0.888161   1.30532   1.57426   1.66667   1.57426   1.30532   0.888161
 0.522393   0.917827  1.28937   1.57426   1.71488   1.65627   1.34723
 0.26094    0.554958  0.917827  1.30532   1.65627   1.88328   1.8543
 0.0981165  0.26094   0.522393  0.888161  1.34723   1.8543    2.28311


```


### Plots
