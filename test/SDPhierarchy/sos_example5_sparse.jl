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

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 1,
                                        issparse = true,
                                        params = Dict(:opt_outlev=>0))

    max_cliques = get_WB5cliques(relax_ctx, problem)
    @assert length(max_cliques) == 2

    logpath = joinpath(testfolder, "$instance")
    ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
    println("Saving file at $logpath")

    primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath,
                                                        save_pbs=true,
                                                        max_cliques=max_cliques);
    @show 1458.8347 dualobj, 1458.8347

    @test primobj ≈ 954.8232 atol=1e-4
    @test dualobj ≈ primobj atol=(mosek_optgap*min(abs(primobj), abs(dualobj)))
end



@testset "WB5 real formulation - order $i - two cliques" for i in 1:2
    instance = "case9"

    problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
    problem = pb_cplx2real(problem_c)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = i,
                                        issparse = true,
                                        params = Dict(:opt_outlev=>0))

    max_cliques = get_case9cliques(relax_ctx, problem)
    @assert length(max_cliques) == 3

    logpath = joinpath(testfolder, "$instance")
    ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
    println("Saving file at $logpath")

    primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath,
                                                        save_pbs=true,
                                                        max_cliques=max_cliques);

    if i == 1
        @show (primobj, dualobj, 1458.60)
        @test primobj ≈ 1458.6003 atol=1e-4
        @test dualobj ≈ primobj atol=(mosek_optgap*min(abs(primobj), abs(dualobj)))
    else
        @show (primobj, dualobj, 1458.60)
        @test primobj ≈ 1458.60 atol=1e-2
        @test dualobj ≈ primobj atol=(mosek_optgap*min(abs(primobj), abs(dualobj)))
    end
end