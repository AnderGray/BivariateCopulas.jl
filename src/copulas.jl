######
# This file is part of the BivariateCopulas.jl package.
#
#   Definition of bivariate copula type, Joint distribution type and methods
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
######

###
#   Notes: This is a module for performing tests with copulas for plotting, sampling and
#   making bivariate joint distributions
#
#   Would it be possible to represent a copula as an expansion of orthogonal
#   polynomials? Could be used as a sentivitivty analysis for dependancy. Something equivilant
#   to disperive Monte Carlo, where not only the correlatin coefficient is changed but the entire
#   dependancy structure. Could this expansion be used in the arithmetic?
#
#   Or can we do it as sums of the Frechet bounds? (Shuffles of M)
#   See: Nelson page 68.
#
#
###

###
#
#   Known bugs:
#                   ->  Differentiating with the interpolators is really ugly
#
###



#import PyPlot: PyPlot.plot

mutable struct copula <: AbstractCopula

    cdf ::  Array{Float64,2}
    density :: Union{Array{Float64,2},Missing}
    func :: Union{Function,Missing}
    param :: Union{Float64,Missing}

    function copula(cdf = missing; density = missing, func = missing, param = missing)

        if (ismissing(cdf) && ismissing(density)) throw(ArgumentError("Cdf and/or density must be provided")) end
        if (ismissing(cdf)); cdf = calcCdf(density); end
        C = new(cdf,density,func,param)
        if (ismissing(density)); density = calcDensity(C); end

        return new(cdf,density, func, param)

    end
end

function (obj::copula)(X,Y, useInterp = false)             # Evaluating the copula

    #if (!issorted(X) || !issorted(Y)) throw(ArgumentError("cdf evaluation request must be given in assending order"));end

    if !useInterp
        if !ismissing(obj.func)
            if !ismissing(obj.param)
                return obj.func(X,Y,obj.param)
            end
            return obj.func(X,Y)
        end
    end

    cdf = interpolate(obj.cdf, BSpline(Linear()))

    uu = range(0,1,length = n)
    sCdf = Interpolations.scale(cdf, uu, uu);
    return sCdf(X,Y);
end

cdf(C:: copula, X, Y) = C(X,Y)

function density(C :: copula, X, Y)

    den = interpolate(C.density, BSpline(Quadratic(Line(OnCell()))))
    uu = range(0,1,length = n)
    sDen = Interpolations.scale(den, uu, uu);
    return sDen(X, Y);
end

function (obj::copula)(X :: ContinuousUnivariateDistribution,Y :: ContinuousUnivariateDistribution)             #Constructing a joint distribution using the copula
    return Joint(X, Y,obj)
end

function sample(C :: AbstractCopula, N ::Int64 = 1; plot = false)

    useInterp = false;      # Set true to use interpolator
    if ismissing(C.func) func = 0; useInterp = true else func = C.func end

    if plot return samplePlot(C,N);end
    if (func == Gau) return CholeskyGaussian(N, C.param) ;end  # Use Cholesky decompostition of the Cov matrix for gaussian copula sampling

    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;

    if func   == perf return hcat(ux,ux); end
    if func   == opp return hcat(ux,1 .- ux); end
    if func   == indep return hcat(ux,y); end

    # if (!ismissing(C.func) && n < 10^5) m = 10^5; end  # If function is provided more acurate sampling, otherwise must use interpolator
    e = sqrt(eps());      ygrid = range(0,stop = 1,length = m);

    #if (C.func == Gau) useInterp = true ;end           # Should you not want to use Cholesky and interpolator has better performance for gaussian

    if func !== 0; conditional = (C(x .+ e/2,ygrid, useInterp) - C(x .- e/2,ygrid, useInterp))./e
    else

    end
    conditional = (C(x .+ e/2,ygrid, useInterp) - C(x .- e/2,ygrid, useInterp))./e
    for i=1:N
        #conditional = (C(x[i] + e/2,ygrid, true) - C(x[i] - e/2,ygrid, true))./e
        #conditional[1] = 0; conditional[end]=1;       # Small errors sometimes happen
        #if (conditional[end] < y[i]) uy[i] = 1;      # Should linearly interpolate between cond[end] and 1
        #elseif (conditional[1] > y[i]) uy[i] = 0;      # Should linearly interpolate between cond[end] and 1
        #else uy[i] = findlast(conditional .<= y[i]) / m; end        # Finds Nothing??
        #yitp_linear = interpolate(conditional, BSpline(Linear()))
        #uy[i] = yitp_linear(y[i] * m);

        #uy[i] = findlast(conditional .<= y[i])/m; end
        if !issorted(conditional[i,:]) conditional[i,:] = sort(conditional[i,:]) end
        inpl = LinearInterpolation(conditional[i,:], ygrid, extrapolation_bc = Flat())
        uy[i] = inpl(y[i])
    end
    return hcat(ux,uy);
end

rand(C :: AbstractCopula, N :: Int64 = 1; plot= false) = sample(C, N, plot=plot)

function CholeskyGaussian(N = 1, correlation = 0)

    Cov = ones(2,2); Cov[1,2] = Cov[2,1] = correlation;
    a = cholesky(Cov).L;
    z = transpose(rand(Normal(),N,2))
    x = a * z;
    u = transpose(cdf.(Normal(),x))

    return hcat(u[:,1],u[:,2])
end

function conditional(C :: AbstractCopula, xVal :: Real; plot = false)   # May also work for joint? xVal will be take from invCdf(M1, u1)

    e = 1/n;        ygrid = range(0,stop=1,length = n);
    if !ismissing(C.func) e = (sqrt(eps())); end

    conditional = (C(xVal + e/2,ygrid) - C(xVal - e/2,ygrid))./e;
    if (plot) plot2(C.cdf, conditional, xVal) end
    return conditional

end

function calcCdf(density :: Array{<:Real,2})
    cdf = zeros(size(density))
    for i = 2:n
        for j = 2:n
            cdf[i,j] =  (density[i,j]/(n^2) + cdf[i, j-1] + cdf[i-1, j] - cdf[i-1, j-1]);
        end
    end
    return cdf
end

function calcDensity(C :: AbstractCopula)

    if ismissing(C.func) return calcDensity(C.cdf); end

    e = sqrt(sqrt(eps()));

    x = range(0 , stop = 1 - e    ,length = n);
    y = range(0 , stop = 1 - e    ,length = n);

    der1 = ( C(x .+ e, y)       - C(x, y)     )/e      # df/dx
    der2 = ( C(x .+ e, y .+ e ) - C(x, y .+ e) )/e     # d(df/dx)/dy
    density = (der2-der1)/e

    return density
end

function calcDensity(cdf :: Array{<:Real,2}, xRange = range(0,1,length = n), yRange = range(0,1,length = n))

    xStep = Float64(xRange.step); yStep = Float64(yRange.step);

    density = zeros(size(cdf));
    for i = 1:n-1
        for j = 1:n-1
            density[i,j] = (cdf[i,j] + cdf[i+1,j+1] - cdf[i+1,j] - cdf[i,j+1])/(xStep*yStep);

            #   -> Possible alternative:
            #der1 = ( cdf[i + 1, j]      - cdf[i, j]     )/(xStep^2)             # df/dx
            #der2 = ( cdf[i + 1, j + 1 ] - cdf[i, j + 1] )/(yStep^2)             # d(df/dx)/dy
            #density[i,j] = (der2 - der1)
        end
    end

    # Problem: How do you get the end values?
    #   Two solutions:

    #   -> If copula, it's symetric (exact answer exept at two endpoints)
    #density[end,:] = reverse(density[:,1]);
    #density[:,end] = reverse(density[1,:]);

    #   -> If joint, it may not be symetric
    density[end,:] = density[end-1,:];
    density[:,end] = density[:,end-1];
    return density
end

function M()
    x = y = range(0,stop=1,length = n);
    return copula(perf(x,y), density = Matrix{Float64}(I, n, n), func = perf);
end

function W()
    x = y = range(0,stop=1,length = n);

    den = zeros(n,n)
    for i = 1:n
        for j = 1:n
            if (i-1) == n-(j-1); den[i,j] = 1; end
        end
    end
    return copula(opp(x,y), density = den , func = opp);
end

function Pi()
    x = y = range(0,stop=1,length = n);
    return copula(indep(x,y),density= ones(n,n), func = indep);
end

function Frank(s = 1)                                     #   Frank copula
    x = y = range(0,stop = 1,length = n);                          #   s>0; 0 for perfect, 1 for indep, inf for oposite
    return copula(F(x,y,s), func = F, param = s);
end

function Clayton(t = 0)                                  #   Clayton copula
    x = y = range(0, stop=1, length = n);                          #   t>-1; -1 for opposite, 0 for indep and inf for perfect
    return copula(Cla(x,y,t), func = Cla, param = t);
end

function Gaussian(corr = 0)
    x = y = range(0,stop=1,length = n);
    cdf = Gau(x,y,corr);
    #cdf[:,end] = cdf[end,:] = x;    # This removes the NaNs from Liams mvNormCdf. We know the marginals must be uniform
    return copula(cdf, func = Gau, param = corr);
end


function τCopula( τ = 0 )   # Imprecise copula from kendal tau

    x = y = range(0,stop = 1,length = n);

    cdf = kenLB(x,y,τ)

    return copula(cdf, func = kenLB, param = τ)

end

KendalCopula(τ = 0) = τCopula(τ)

function ρCopula( ρ = 0 ) # Imprecise copula from Spearman rho
    
    x = y = range(0,stop = 1,length = n);

    cdf = spearLB(x,y,ρ)

    return copula(cdf, func = spearLB, param = ρ)
end

SpearmanCopula(ρ = 0 ) = ρCopula( ρ )

indep(X, Y) = [x*y for x in X, y in Y];
perf(X, Y)  = [min(x,y) for x in X, y in Y];
opp(X, Y)   = [max(x+y-1,0) for x in X, y in Y];
F(X,Y,s = 1)      = [log(1+(s^x-1)*(s^y-1)/(s-1))/log(s) for x in X, y in Y]
Cla(X,Y, t = 0) = [max((x^(-t)+y^(-t)-1)^(-1/t),0) for x in X, y in Y]
Gau(X,Y, corr = 0)  = [bivariate_cdf(quantile.(Normal(),x),quantile.(Normal(),y), corr) for x in X, y in Y];


#kenLB(x, y, τ) = max(0, x+y-1, 0.5*( (x+y)-sqrt( (x-y)^2)+1-τ ) )
kenUB(x, y, τ) = min(x, y, 0.5*( (x+y-1) + sqrt( (x+y-1)^2 +1+τ )) )
kenLB(x, y, τ) = [x - kenUB(x, 1 - y, - τ) for x in x, y in y]

#spearϕ(a, b) = 1/6 * ( (max(0, 9*b + 3*sqrt( 9*b^2 - 3*a^6)))^(1/3) + ( max(0,9*b - 3*sqrt(9*b^2 - 3*a^6)))^(1/3) )

function spearϕ(a, b)
    A = 9*b 
    B = max(9*b^2 - 3*a^6, 0) 
    C = (max(0, A + 3*sqrt(B)))^(1/3)
    D = (max(0, A - 3*sqrt(B)))^(1/3)
    return 1/6 * (C + D)
end

#spearLB(x, y, ρ) = max(0, x + y - 1, (x + y)/2 - spearϕ(x - y, 1 - ρ))
spearUB(x, y, ρ) = max( opp(x,y), min(x, y, (x + y -1)/2 + spearϕ(x + y - 1, 1 + ρ)))
spearLB(x, y, τ) = [x - spearUB(x, 1 - y, - τ) for x in x, y in y]


function GauCopulaDen(corr = 0, X = range(bOt, tOp, length=n), Y = range(bOt, tOp, length=n))

    COR = [1. corr; corr 1.];
    I = [1. 0; 0 1.]
    den = zeros(length(X),length(Y))
    for i = 1:length(X)
        for j = 1:length(Y)
            Vec = [invCdf(Normal(),X[i]); invCdf(Normal(),Y[j])]
            den[i,j] = 1/sqrt(det(COR)) * exp(-1/2 * Vec' * (inv(COR) - I) *Vec)
        end
    end
    return den
end

Arch1(X,Y,theta = 0) = [1 - ((1-x)^theta + (1-y)^theta - (1-x)^theta * (1-y)^theta )^theta for x in X, y in Y];

mvNormCdf(X,Y,cor) = [bivariate_cdf(x,y,cor) for x in X, y in Y];

mutable struct Joint <: AbstractJoint

    marginal1 :: ContinuousUnivariateDistribution
    marginal2 :: ContinuousUnivariateDistribution
    copula  :: AbstractCopula

    xRange
    yRange

    cdf ::  Union{Array{Float64,2},Missing}     # Or should the density/ cdf be calculated on the fly from the copula?
    density :: Union{Array{Float64,2},Missing}

    function Joint(marginal1, marginal2, copula)

        xRange = range(invCdf(marginal1,bOt), stop=invCdf(marginal1,tOp), length = n)
        yRange = range(invCdf(marginal2,bOt), stop=invCdf(marginal2,tOp), length = n)

        j = new(marginal1, marginal2, copula, xRange, yRange, missing)
        cdf = j(xRange,yRange)
        j = new(marginal1, marginal2, copula, xRange, yRange, cdf)
        density = calcDensity(j)                                    # Looks really ugly when using the interpolator for a marginal

        return new(marginal1, marginal2, copula, xRange, yRange, cdf, density)

    end
end

function (obj :: Joint)(X,Y)    # Evaluating the joint cdf
                                # BiLinear interpolator if cdf is provided

    if !ismissing(obj.copula.func)
        return obj.copula(cdf.(obj.marginal1,X), cdf.(obj.marginal2,Y))  # This is the case where we have the functions
    end
    if ismissing(obj.cdf)
        return obj.copula(cdf.(obj.marginal1,X), cdf.(obj.marginal2,Y))
    end

    cdf1 = interpolate(obj.cdf, BSpline(Quadratic(Line(OnCell()))))
    sCdf1 = Interpolations.scale(cdf1, obj.xRange, obj.yRange);

    return sCdf1(X, Y);

    # Should use new scaling from interpolatorrs
end

cdf(J::Joint, X, Y) = J(X, Y)

function conditional(J :: Joint, xVal :: Real; plot = false)   # May also work for joint? xVal will be take from invCdf(M1, u1)

    e = 1/n;        ygrid = J.yRange;
    if !ismissing(J.copula.func) e = (sqrt(eps())); end

    conditional = (J(xVal + e/2,ygrid) - J(xVal - e/2,ygrid))./e;

    conditional = conditional/pdf(J.marginal1,xVal)

    if (plot) plot2(J.cdf, conditional, xVal) end
    return conditional
    
end

function density2(J :: Joint, X = J.xRange, Y = J.yRange)

    den = interpolate(J.density, BSpline(Quadratic(Line(OnCell()))))
    sDen = Interpolations.scale(den, X, Y);
    return sDen(X,Y);

end

function density(J :: Joint,X = J.xRange, Y = J.yRange)        # An alternative way to construct the Joint. fxy / fx*fy = c[Fx,Fy]. Where c is the copula Density
                                                                # More error found against Multivariate gaussian pdf. Due to interpolation error in the copula density.
    denCop = density(J.copula, cdf.(J.marginal1,X), cdf.(J.marginal2,Y))
    return denCop .* [pdf(J.marginal1, x) * pdf(J.marginal2, y) for x in X, y in Y]

end

function sample(J :: AbstractJoint, N :: Int64 =1)

    copulaSamples = sample(J.copula,N)

    x = invCdf.(J.marginal1,copulaSamples[:,1]);
    y = invCdf.(J.marginal2,copulaSamples[:,2]);
    return hcat(x,y)

end

rand(J :: AbstractJoint, N :: Int64 = 1) = sample(J :: AbstractJoint, N)

function invCdf(X :: Sampleable{Univariate}, u)
    return quantile.(X,u)                           # We will also need the interpolartor here
end

function calcDensity(J :: AbstractJoint)        

    #return calcDensity(J.cdf, J.xRange, J.yRange)

    if ismissing(J.copula.func) return calcDensity(J.cdf,J.xRange,J.yRange) end

    e = sqrt(sqrt(eps()));

    #x = J.xRange - e;f
    #y = J.yRange - e;

    x = range(J.xRange[1], stop = J.xRange[end]-e, length = length(J.xRange))
    y = range(J.xRange[1], stop = J.yRange[end]-e, length = length(J.yRange))

    der1 = ( J(x .+ e, y)       - J(x, y)     )/e     # df/dx
    der2 = ( J(x .+ e, y .+ e ) - J(x, y .+ e) )/e     # d(df/dx)/dy
    density = (der2-der1)/e                            #d^2f/dxdy

    return density

end

M(M1::ContinuousUnivariateDistribution, M2::ContinuousUnivariateDistribution) = M()(M1,M2)
W(M1::ContinuousUnivariateDistribution, M2::ContinuousUnivariateDistribution) = W()(M1,M2)
Pi(M1::ContinuousUnivariateDistribution, M2::ContinuousUnivariateDistribution) = Pi()(M1,M2)
Frank(M1::ContinuousUnivariateDistribution, M2::ContinuousUnivariateDistribution, s = 1) = Frank(s)(M1,M2)
Clayton(M1::ContinuousUnivariateDistribution, M2::ContinuousUnivariateDistribution, t = 0) = Clayton(t)(M1,M2)
Gaussian(M1::ContinuousUnivariateDistribution, M2::ContinuousUnivariateDistribution, corr = 0) = Gaussian(corr)(M1,M2)

struct marginal <: AbstractMarginal

    cdf :: Array{Float64,1}
    range
    cdfInp
    invCdfInp

    # Could have density and a calcDensity function for pdf

    function marginal(cdf :: Array{Float64,1}, range)

        if (length(cdf) !== length(range)) throw(ArgumentError("Cdf and range must be the same length")) end

        inpl = interpolate(cdf, BSpline(Quadratic(Line(OnCell()))));
        cdfInp = Interpolations.scale(inpl, range);

        # Error here and also in differention
        invCdfInp = LinearInterpolation(cdf, range, extrapolation_bc = Flat())

        new(cdf, range, cdfInp, invCdfInp)
    end
end

function (obj:: marginal)(x)
    return cdf(obj, x)
end

function cdf(D :: AbstractMarginal, x :: Real)

    if x < D.range[1] return 0; end
    if x > D.range[end] return 1; end

    return D.cdfInp(x)

end


function invCdf(D :: AbstractMarginal, u :: Real)

    return D.invCdfInp(u)

end



quantile(D :: marginal, u) = invCdf(D, u);

function sample(D :: marginal, N = 1)
    return invCdf(D, rand(N))
end

rand(D::marginal, N=1) = sample(D, N)


function Base.show(io::IO, z::copula)

    statement1 = "Arbitrary"
    statement2 = ""

    if (!ismissing(z.func))
        func = z.func
        if (func == indep); func ="π";end
        if (func == perf); func = "M";end
        if (func == opp); func = "W";end
        statement1 = "$(func)"
    end

    if (!ismissing(z.param))
        func = z.func
        parName = "par"
        if (func == Gau) parName = "r";end
        if (func == F) parName = "s";end
        if (func == Cla) parName = "t";end
        statement2 = "$parName=$(z.param)"
    end

    print(io, "Copula ~ $statement1($statement2)");
end

function Base.show(io::IO, z::Joint)

    statement1 = "Arbitrary"
    statement2 = ""

    if (!ismissing(z.copula.func))
        func = z.copula.func
        if (func == indep); func ="π";end
        if (func == perf); func = "M";end
        if (func == opp); func = "W";end
        statement1 = "$(func)"
    end

    if (!ismissing(z.copula.param))
        func = z.copula.func
        parName = "par"
        if (func == Gau) parName = "r";end
        if (func == F) parName = "s";end
        if (func == Cla) parName = "t";end
        if (func == spearLB) parName = "ρ";end 
        if (func == kenLB) parName = "τ";end
        statement2 = "$parName=$(z.copula.param), "
    end

    print(io, "Joint ~ $statement1( $statement2$(z.marginal1), $(z.marginal2) )");
end


function Base.show(io::IO, z::marginal)

    print(io, "Distribution ~ ( range=[$(z.range[1]),$(z.range[end])])");
end


function linInterp(A, x)

    ind0 = findlast(A .<= x);
    x0 = A[ind0];   x1 = A[ind0+1]
    y0 = ind0/n;     y1 = (ind0+1)/n

    return y0 + (x-x0) * (y1-y0)/(x1-x0)
end

