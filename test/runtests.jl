using MathProgComplex, DataStructures, Memento, Printf, Pkg
MPC = MathProgComplex

using DataStructures, Memento, Test

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# suppress all messages less important than errors
setlevel!(getlogger(MathProgComplex), "error")

# @testset "MathProgComplex" begin
    include("basetypes.jl")
    include("cplx2real.jl")
    include("dat_importexport.jl")
    include("polynomial_types.jl")
    include("jump_export.jl")
# end
