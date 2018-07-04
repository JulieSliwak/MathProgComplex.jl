using MathProgComplex, DataStructures

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
include("basetypes.jl")
include("cplx2real.jl")
include("dat_importexport.jl")
include("polynomial_types.jl")
include("jump_export.jl")
