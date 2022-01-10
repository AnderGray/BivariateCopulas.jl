######
# This file is part of the BivariateCopulas.jl package.
#
#   Module definition
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Authors: Ander Gray & Liam T. Henry
#                                           Email:  ander.gray@liverpool.ac.uk
######

module BivariateCopulas

using Reexport
using Requires

@reexport using Distributions             #   https://github.com/JuliaStats/Distributions.jl
using LinearAlgebra                       #   Main Library
using Interpolations                      #   https://github.com/JuliaMath/Interpolations.jl

import Distributions: cdf, quantile, rand

import Base: rand

global n = 200                       # Number of descretizations
global pn = 50                       # Number of descretization for plotting
global bOt = 0.0001                    # smallest quantile to use if left tail is unbounded. Should we use these to define the copula ranges?
global tOp = 0.9999                    # largest quantile to use if right tail is unbounded
#global bOt = 0;
#global tOp = 1;

abstract type AbstractCopula <: Real end
abstract type AbstractJoint <: Real end
abstract type AbstractMarginal <: ContinuousUnivariateDistribution end

include("archimedean/archimedean.jl")
include("copulas.jl")
include("NormalDistribution.jl")

function __init__()
    @require PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        @require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" begin
            using PyCall
            using PyPlot
            import PyPlot: plot
            include("plotting.jl")
            return using3D()
        end
    end
end

export copula,
    density,
    cdf,
    sample,
    rand,
    CholeskyGaussian,
    conditional,
    calcCdf,
    calcDensity,
    M,
    W,
    Pi,
    Frank,
    Clayton,
    Gaussian,
    norm_cdf,
    bivariate_cdf,
    mvNormCdf,
    Joint,
    invCdf,
    marginal,
    quantile,
    ρCopula,
    SpearmanCopula,
    KendalCopula,
    τCopula,
    samplePlot,
    plotDensity,
    plotDen,
    plotCdf,
    plot,
    plotContourCdf,
    plotContourDen,
    scatter,
    rotate90,
    rotate180,
    rotate270,
    Clayton,
    Independence,
    rosenblatt,
    inverse_rosenblatt,
    τ

end #BivariateCopulas
