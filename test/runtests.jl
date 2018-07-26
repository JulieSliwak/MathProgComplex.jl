using MathProgComplex, DataStructures

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# # Polynomial optimization problems
# @testset "Polynomial optim" begin
#     include("basetypes.jl")
#     include("cplx2real.jl")
#     include("dat_importexport.jl")
#     include("polynomial_types.jl")
#     include("jump_export.jl")
# end

const mosek_optgap = 1.0e-6

struct RelaxSol
    instance::String
    relaxation_order::Int
    primal_solvalue::Float64
    rel_opt_gap::Float64
end

const OPFsols = OrderedDict(("WB2", 1)             => RelaxSol("WB2",             1, 885.7145992, 1.00E-03),
                            ("LMBM3", 1)           => RelaxSol("LMBM3",           1, 386.4164287, 1.00E-03),
                            ("WB5", 1)             => RelaxSol("WB5",             1, 954.8231857, 1.00E-03),
                            ("case6ww", 1)         => RelaxSol("case6ww",         1, 2986.033688, 1.00E-03),
                            ("case9mod", 1)        => RelaxSol("case9mod",        1, 1319.616781, 1.00E-03),
                            ("case9", 1)           => RelaxSol("case9",           1, 1458.83467, 1.00E-03),
                            ("case14", 1)          => RelaxSol("case14",          1, 5371.499451, 1.00E-03),
                            ("case24_ieee_rts", 1) => RelaxSol("case24_ieee_rts", 1, 50711.39264, 1.00E-03),
                            ("case_ieee30", 1)     => RelaxSol("case_ieee30",     1, 5927.58577, 1.00E-03),
                            ("case30pwl", 1)       => RelaxSol("case30pwl",       1, 71.99999998, 1.00E-03),
                            ("case30", 1)          => RelaxSol("case30",          1, 316.48926, 1.00E-03),
                            ("WB2", 2)             => RelaxSol("WB2",             2, 890.07615, 1.00E-03),
                            ("LMBM3", 2)           => RelaxSol("LMBM3",           2, 386.4164767, 1.00E-03),
                            ("WB5", 2)             => RelaxSol("WB5",             2, 1146.478816, 1.00E-03),
                            ("case6ww", 2)         => RelaxSol("case6ww",         2, 2986.034849, 1.00E-03),
                            ("case9mod", 2)        => RelaxSol("case9mod",        2, 1320.558097, 1.00E-03),
                            ("case9", 2)           => RelaxSol("case9",           2, 1458.848541, 1.00E-03),
                            ("case14", 2)          => RelaxSol("case14",          2, 5371.500624, 1.00E-03))

# SDP hierarchy problems
@testset "SDP hierarchy" begin

    include(joinpath("SDPhierarchy", "sos_example1.jl"))
    include(joinpath("SDPhierarchy", "sos_example2.jl"))
    include(joinpath("SDPhierarchy", "sos_example5_sparse.jl"))
    include(joinpath("SDPhierarchy", "sos_example6_matpower_rankrel.jl"))

    # include(joinpath("SDPhierarchy", "sos_example3.1_dense_WB2.jl"))
    # include(joinpath("SDPhierarchy", "sos_example3.2_dense_WB2.jl"))
    # include(joinpath("SDPhierarchy", "sos_example4.1_dense_WB5.jl"))
    # include(joinpath("SDPhierarchy", "sos_example4.2_dense_WB5.jl"))
    # include(joinpath("SDPhierarchy", "sos_example7_matpower_ordre2.jl"))
end