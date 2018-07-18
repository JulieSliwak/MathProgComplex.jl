using Base.Test, MathProgComplex, OPFInstances

## Beware, rank relaxation value is for problem without objective constant term

@testset "WB2 real formulation - order 1" begin
    sols = OrderedDict( ["WB2"            => (2, 885.71,   905.73, true),
                        # "WB3"            => (1, 417.25,   417.25, false),
                        "LMBM3"          => (1, 386.42,   386.42, false),
                        "WB5"            => (2, 954.82,   1146.4, true),
                        "case6ww"        => (1, 2986,     2986,   false),
                        "case9"          => (2, 373.8,    1458.8, true),
                        "case9mod"       => (2, 234.6,    1320.4, true),
                        "case14"         => (2, 721.5,    5371.5, true),
                        # "case22loop"     => (1, 4538.8,   4538.8, false), ## Absent from data repo...
                        "case30"         => (2, 268.915,  316.49, true)]) #,
                        # "case39"         => (1, 1887.2,   1887.2, false)) #,
                        # "case39mod1"     => (2, 773.36,   942.34, true),
                        # "case39mod2"     => (1, 940.34,   940.34, false),
                        # "case57"         => (1, 25338,    25338,  false),
                        # "case89pegase"   => (1, 5817.6,   5817.6, false)) #,
                        # "case118"        => (1, 86298,    86298,  false),
                        # "case118mod"     => (1, 86079,    86079,  false),
                        # "case300"        => (2, 254840.4, 475430, true),
                        # "case300mod"     => (1, 287950,   287950, false),
                        # "case1354pegase" => (2, 74053,    74053,  true))

    testfolder = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "tests", "sos_example6")

    @testset "instance $instance, (rmeqs=$rmeqs)" for (instance, (dcv, obj_rankrel, obj_opt, lackconstant)) in sols, rmeqs in Set([false])

        mosek_optgap = 1e-4

        info("--> Working on $instance")
        # OPFpbs = load_OPFproblems(MatpowerInput, joinpath("..", "data", "data_Matpower", "matpower", instance*".m"))
        # problem_c = build_globalpb!(OPFpbs)

        problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
        cstobj = 0
        (lackconstant) && haskey(problem_c.objective, Exponent()) && (cstobj = problem_c.objective[Exponent()])

        problem = pb_cplx2real(problem_c)

        logpath = joinpath(testfolder, "$instance")
        ispath(logpath) && rm(logpath, recursive=true); mkpath(logpath)
        println("Saving file at $logpath")

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = 1,
                                            params = Dict(:opt_outlev=>0))

        primobj, dualobj = run_hierarchy(problem, relax_ctx, logpath, save_pbs=true);
        @show (primobj, dualobj, obj_rankrel, obj_opt, cstobj)

        dualgaptol::Float64 = mosek_optgap*min(abs(primobj), abs(dualobj))
        objtol::Float64 = 1e-4*abs(obj_rankrel + cstobj)

        objandctr::Float64 = obj_rankrel + cstobj

        @test primobj ≈ objandctr atol=objtol
        @test dualobj ≈ primobj atol=dualgaptol
    end
end
