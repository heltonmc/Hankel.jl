module Hankel
import Bessels
import SpecialFunctions
using GSL: GSL
using Roots: Roots
import LinearAlgebra: mul!, ldiv!, dot
import Base: *, \
using ChainRulesCore

export QDHT, integrateK, integrateR, onaxis, symmetric, Rsymmetric

include("utils.jl")
include("qdht.jl")
include("chainrules.jl")

besselj(v::Real, x::Real) = besselj(promote(v, x)...)
besselj(v::T, x::T) where T <: Union{Float16, Float32, Float64, Int32, Int64} = Bessels.besselj(v, x)
besselj(v::T, x::T) where T <: AbstractFloat = SpecialFunctions.besselj(v, x)

gamma(x::T) where T <: Union{Float16, Float32, Float64, Int32, Int64} = Bessels.gamma(x)
gamma(x::T) where T <: AbstractFloat = SpecialFunctions.gamma(x)

const J₀₀ = besselj(0, 0)

end
