using Base.Test, MathProgComplex


## Application of the Moment–SOS Approach to Global Optimization of the OPF Problem
# C. Josz, J. Maeght, P. Panciatici, and J. Ch. Gilbert
# https://arxiv.org/pdf/1311.6370.pdf
# Data from table 3, p. 13.

@testset "WB5 real formulation - dense - order 2" begin
    sols = SortedDict(-30.80 => (2, 945.83, 945.83),
                      -20.51 => (2, 1146.48, 954.82),
                      -10.22 => (2, 1209.11, 963.83),
                       00.07 => (2, 1267.79, 972.85),
                       10.36 => (2, 1323.86, 981.89),
                       20.65 => (2, 1377.97, 990.95),
                       30.94 => (2, 1430.54, 1005.13),
                       41.23 => (2, 1481.81, 1033.07),
                       51.52 => (2, 1531.97, 1070.39),
                       61.81 => (1, 1114.90, 1114.90))

    testfolder = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "tests")
    @testset "q5min=$q5min, rmeqs=$rmeqs" for (q5min, (dmax, obj_fixedphase, obj_rankrel)) in sols, rmeqs in Set([true])

        # info("Working on WB5, no eqs, q5min=$q5min, d=$d")
        problem = buildPOP_WB5(q5min=q5min, rmeqs=rmeqs)

        logpath = joinpath(testfolder, "WB5_q5min_$(q5min)_ordercv")
        ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
        println("Saving file at $logpath")

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = dmax)

        primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath, save_pbs=true);
        @show primobj, dualobj, obj_fixedphase
        @test primobj ≈ obj_fixedphase atol=1e-2
        @test dualobj ≈ primobj atol=mosek_optgap*min(abs(primobj), abs(dualobj))
    end
end
