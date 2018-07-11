using Base.Test
using MathProgComplex

@testset "WB2 real formulation, fixed phase, order 2" begin
    sols = SortedDict(0.976 => (2, 905.76, 905.76),
                    0.983 => (2, 905.73, 903.12),
                    0.989 => (2, 905.73, 900.84),
                    0.996 => (2, 905.73, 898.17),
                    1.002 => (2, 905.73, 895.86),
                    1.009 => (2, 905.73, 893.16),
                    1.015 => (2, 905.73, 890.82),
                    1.022 => (3, 905.73, 888.08),
                    1.028 => (3, 905.73, 885.71),
                    1.035 => (2, 882.97, 882.97))

    testfolder = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "tests")
    @testset "v2max=$v2max, rmeqs=$rmeqs, addball=$addball" for (v2max, (dmax, obj_fixedphase, obj_rankrel)) in sols, rmeqs in Set([true, false]), addball in Set([true, false])

        # info("Working on WB2, no eqs, v2max=$v2max, d=$d")
        problem = buildPOP_WB2(v2max=v2max, rmeqs=rmeqs, setnetworkphase=true, addball=addball)

        logpath = joinpath(testfolder, "WB2_v2max_$(v2max)_rankrel")
        ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
        println("Saving file at $logpath")

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = dmax)

        primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath, save_pbs=true);
        @show primobj, dualobj, obj_fixedphase
        if (0.983 ≤ v2max ≤ 1.028)
            @test_broken primobj ≈ obj_fixedphase atol=1e-2
            @test_broken dualobj ≈ primobj atol=mosek_optgap*min(abs(primobj), abs(dualobj))
        else
            @test primobj ≈ obj_fixedphase atol=1e-2
            @test dualobj ≈ primobj atol=mosek_optgap*min(abs(primobj), abs(dualobj))
        end
    end
end