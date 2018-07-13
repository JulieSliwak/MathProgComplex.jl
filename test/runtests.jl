using MathProgComplex, DataStructures

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# # Polynomial optimization problems
include("basetypes.jl")
include("cplx2real.jl")
include("dat_importexport.jl")
include("polynomial_types.jl")
include("jump_export.jl")

const mosek_optgap = 1e-6

# SDP hierarchy problems
# include(joinpath("SDPhierarchy", "sos_example1.jl"))
# include(joinpath("SDPhierarchy", "sos_example2.jl"))
# include(joinpath("SDPhierarchy", "sos_example5_sparse.jl"))
# include(joinpath("SDPhierarchy", "sos_example6_matpower_rankrel.jl"))

# include(joinpath("SDPhierarchy", "sos_example3.1_dense_WB2.jl"))
# include(joinpath("SDPhierarchy", "sos_example3.2_dense_WB2.jl"))
# include(joinpath("SDPhierarchy", "sos_example4.1_dense_WB5.jl"))
# include(joinpath("SDPhierarchy", "sos_example4.2_dense_WB5.jl"))
# include(joinpath("SDPhierarchy", "sos_example7_matpower_ordre2.jl"))