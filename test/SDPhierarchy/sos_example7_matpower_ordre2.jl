using Base.Test, MathProgComplex, OPFInstances

## Beware, rank relaxation value is for problem without objective constant term

@testset "WB2 real formulation - order 2" begin
    sols = OrderedDict("WB2"        => (2, 885.71, 905.73, true),
                    #    "WB3"        => (1, 417.25, 417.25, false),
                       "LMBM3"      => (1, 386.42, 386.42, false),
                       "WB5"        => (2, 954.82, 1146.4, true),
                       "case6ww"    => (1, 2986, 2986, false),
                       "case9"      => (2, 373.8, 1458.8, true),
                       "case9mod"   => (2, 234.6, 1320.4, true),
                       "case14"     => (2, 721.5, 5371.5, true))
                    # #    "case22loop" => (1, 4538.8, 4538.8, false), ## Absent in data repo...
                    #    "case30"     => (2, 268.915, 316.49, true))

    for (key, val) in sols
        (val[1] == 1) && delete!(sols, key)
    end

    testfolder = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "tests", "sos_example6")

    @testset "instance $instance, (rmeqs=$rmeqs)" for (instance, (dcv, obj_rankrel, obj_opt, lackconstant)) in sols, rmeqs in Set([false])

        info("--> Working on $instance")
        # OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", instance*".m"))
        # problem_c = build_globalpb!(OPFpbs)

        problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
        cstobj = 0
        # (lackconstant) && haskey(problem_c.objective, Exponent()) && (cstobj = problem_c.objective[Exponent()])

        problem = pb_cplx2real(problem_c)

        logpath = joinpath(testfolder, "$instance")
        ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
        println("Saving file at $logpath")

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            symmetries=[PhaseInvariance],
                                            d = 2)

        primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath, save_pbs=true);
        @show (primobj, dualobj, obj_rankrel, obj_opt, cstobj)

        @test primobj ≈ obj_opt + cstobj atol=1e-1
        @test dualobj ≈ primobj atol=mosek_optgap*min(abs(primobj), abs(dualobj))
    end
end
