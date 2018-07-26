using Base.Test, MathProgComplex, OPFInstances

## Application of the Moment–SOS Approach to Global Optimization of the OPF Problem
# C. Josz, J. Maeght, P. Panciatici, and J. Ch. Gilbert
# https://arxiv.org/pdf/1311.6370.pdf
# Data from table 3, p. 13.

testfolder = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "tests", "sos_example5")

@testset "WB5 real formulation - order 1 - two cliques" begin
    instance = "WB5"

    problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
    problem = pb_cplx2real(problem_c)

    logpath = joinpath(testfolder, "$instance")

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 1,
                                        issparse = true,
                                        params = Dict(:opt_outlev=>0,
                                                      :opt_logpath=>logpath))

    max_cliques = get_WB5cliques(relax_ctx, problem)
    @assert length(max_cliques) == 2

    primobj, dualobj = run_hierarchy(problem, relax_ctx, save_pbs=true,
                                                         max_cliques=max_cliques);


    ε_rel = OPFsols[("WB5", 1)].rel_opt_gap
    ε_abs = max(primobj, dualobj) * ε_rel

    @show primobj, dualobj,  OPFsols[("WB5", 1)].primal_solvalue, ε_abs

    @test primobj ≈ OPFsols[("WB5", 1)].primal_solvalue atol=ε_abs
    @test dualobj ≈ primobj atol=ε_abs
end



@testset "case9 real formulation - order $i - three cliques" for i in 1:2
    instance = "case9"

    problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
    problem = pb_cplx2real(problem_c)

    logpath = joinpath(testfolder, "$instance")

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = i,
                                        issparse = true,
                                        params = Dict(:opt_outlev=>0,
                                                      :opt_logpath=>logpath))

    max_cliques = get_case9cliques(relax_ctx, problem)
    @assert length(max_cliques) == 3

    primobj, dualobj = run_hierarchy(problem, relax_ctx, save_pbs=true,
                                                         max_cliques=max_cliques);


    ε_rel = OPFsols[("case9", i)].rel_opt_gap
    ε_abs = max(primobj, dualobj) * ε_rel

    @show primobj, dualobj,  OPFsols[("case9", i)].primal_solvalue, ε_abs

    @test primobj ≈ OPFsols[("case9", i)].primal_solvalue atol=ε_abs
    @test dualobj ≈ primobj atol=ε_abs
end