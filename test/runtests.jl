using MathProgComplex, DataStructures

include(joinpath("test_pbs", "WB2real.jl"))

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
include("basetypes.jl")
include("cplx2real.jl")
# include("dat_export.jl")
# include("polynomial_types.jl")
# include("jump_export.jl")
