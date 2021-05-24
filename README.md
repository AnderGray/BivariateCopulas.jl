# 2-copulas julia

An implementation of bivariate copulas and bivariate distributions in Julia.

This module has the following capabilities:
* Construction of continous bivariate distributions with Distributions.jl
* Evaluation of cdf and density
* Conditioning (cdfs)
* Sampling
* Bivariate scatter plots, cdf and density plots and contour plots
* Common copula families:
  * π (independence)
  * M (maximal correlation)
  * W (maximal negative correlation)
  * Gaussian (correlation coefficient r; r = -1 for W, r = 0 for π, r = 1 for M)
  * Frank (s is real; -inf for W, 0 for π, inf for M)
  * Clayton (t>=-1; -1 for W, 0 for π and inf for M)

Any possible multivariate dependence can be encoded in a [copula](https://en.wikipedia.org/wiki/Copula_(probability_theory)). Copulas are joint cdfs with standard uniform marginals, and are a way to model dependency independently of marginal distributions. Most probabilistic dependence problems can be reduced to an analysis of copulas.

### Authors

- Ander Gray, Institute for Risk and Uncertainty, University of Liverpool
- Liam T. Henry, Clarus Financial Technology, Belfast

## Installation

This is a registered julia package:

```Julia
julia> ]
pkg> add BivariateCopulas
```

or the most up to date version:

```Julia
julia> ]
pkg> add https://github.com/AnderGray/BivariateCopulas.jl#master
```

## Sklar's theorem

A bivariate copula (2-copula) C is a function C:[0,1]<sup>2</sup> -> [0,1] with the following properties: 

<img src="https://imgur.com/8mDL4fr.png" data-canonical-src="https://imgur.com/8mDL4fr.png" width="500" />

Three important 2-copulas are:

<img src="https://imgur.com/bCIWIDX.png" data-canonical-src="https://imgur.com/bCIWIDX.png" width="400" />

with W and M being bounds on all 2-copulas: W ≤ C ≤ M

Copulas are mainly used in dependence modelling, and can be used to construct any continous multivariate distribution function (df) given their univariate marginals. This is enabled by a theorem from Sklar:


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

### Sklar's theorem

Use any two continous distribution from Distributions.jl may be used to make a bivariate distribution with C

```Julia
julia> C1 = Gaussian(0.8);

julia> j1 = C1(Normal(0, 1), Normal(0, 1)) # Correlated normal
Joint ~ Gau( r=0.8, Normal{Float64}(μ=0.0, σ=1.0), Normal{Float64}(μ=0.0, σ=1.0) )

julia> C2 = M();

julia> j2 = C2(Beta(2, 4), Beta(4, 2)) # Maximally correlated betas
Joint ~ M( Beta{Float64}(α=2.0, β=4.0), Beta{Float64}(α=4.0, β=2.0) )

julia> j3 = W(Normal(-1, 20), Beta(4, 2)) # Minimally correlated normal and beta
Joint ~ W( Normal{Float64}(μ=-1.0, σ=20.0), Beta{Float64}(α=4.0, β=2.0) )

```
### CDF and Density for C and Joints

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

julia> j1 = cop(Beta(2, 4), Normal(5, 3.4));

julia> j1(0.1:0.2:0.8, 5.1:0.2:5.8) # cdf
4×4 Array{Float64,2}:
 0.0806019  0.0807922  0.0809446  0.0810654
 0.388167   0.398055   0.407184   0.415547
 0.5023     0.523831   0.544962   0.565593
 0.511647   0.535043   0.558309   0.581364

julia> density(j1, 0.1:0.2:0.8, 5.1:0.2:5.8) # cdf
4×4 Array{Float64,2}:
 0.0461003   0.0380492   0.0311039   0.0251831
 0.398346    0.390902    0.379929    0.365733
 0.128563    0.142985    0.157505    0.171839
 0.00374385  0.00473385  0.00592841  0.00735337

```

### Sampling

``` Julia
julia> rand(cop, 10^4) # sample a copula
10000×2 Array{Float64,2}:
 0.342484   0.293974
 0.970305   0.892417
 0.344887   0.670971
 0.177868   0.091284
 ⋮          
 0.577286   0.439561
 0.560594   0.405275
 0.259065   0.0325417
 0.424095   0.389324

julia> rand(j1, 10^4) # sample a bivariate distribution
10000×2 Array{Float64,2}:
 0.0737024  -0.791622
 0.330924    4.26197
 0.14506     5.36363
 0.353489    7.28413
 0.530191   10.199
 ⋮          
 0.366434    5.22849
 0.109036    1.30304
 0.555351    6.06244
 0.653661   10.6358
 0.242393    4.85661

```
### Plots

#### Scatter plots
``` Julia
julia> samps = rand(cop, 10^4);

julia> scatter(samps)
```

<img src="https://imgur.com/eot44fc.png" data-canonical-src="https://imgur.com/eot44fc.png" width="700" />

``` Julia
julia> samps2 = rand(j1, 10^4);

julia> scatter(samps2)
```

<img src="https://imgur.com/0cTlomm.png" data-canonical-src="https://imgur.com/0cTlomm.png" width="700" />


#### CDF plots
``` Julia
julia> using3D()
julia> plot(cop)
``` 
<img src="https://imgur.com/2MioWzX.png" data-canonical-src="https://imgur.com/2MioWzX.png" width="700" />

``` Julia
julia> plot(j1)
```
<img src="https://imgur.com/zJ4NSWg.png" data-canonical-src="https://imgur.com/zJ4NSWg.png" width="700" />

#### Density plots
``` Julia
julia> c1 = Frank(0.8);

julia> plotDen(c1)
``` 
<img src="https://imgur.com/Tf6G0cr.png" data-canonical-src="https://imgur.com/Tf6G0cr.png" width="700" />


``` Julia

julia> j1 = c1(Beta(2, 5), Beta(2, 5));

julia> plotDen(j1)
```
<img src="https://imgur.com/VKPSZZ1.png" data-canonical-src="https://imgur.com/VKPSZZ1.png" width="700" />

#### Contour plots

``` Julia
julia> plotContourCdf(j1)
```

<img src="https://imgur.com/AHUxsT5.png" data-canonical-src="https://imgur.com/AHUxsT5.png" width="700" />

``` Julia
julia> plotContourDen(j1)
```

<img src="https://imgur.com/OxDlxQq.png" data-canonical-src="https://imgur.com/OxDlxQq.png" width="700" />


Bibliography
---

* [*Ferson, S., R. Nelsen, J. Hajagos, D. Berleant, J. Zhang, W.T. Tucker, L. Ginzburg and W.L. Oberkampf. 2004. Dependence in Probabilistic Modeling, Dempster-Shafer Theory, and Probability Bounds Analysis. Sandia National Laboratories, SAND2004-3072, Albuquerque, NM*](https://www.osti.gov/servlets/purl/1427286)

* *Nelsen, Roger B. An introduction to copulas. Springer Science & Business Media, 2007*

* *Joe, Harry. Multivariate models and multivariate dependence concepts. CRC Press, 1997*

