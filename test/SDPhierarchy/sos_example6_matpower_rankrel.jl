using Base.Test, MathProgComplex, OPFInstances

## Beware, rank relaxation value is for problem without objective constant term

@testset "ACOPFs - real formulation - order 1" begin
    # sols = OrderedDict( ["WB2"            => (2, 885.71,   905.73, true),
    #                     # "WB3"            => (1, 417.25,   417.25, false),
    #                     "LMBM3"          => (1, 386.42,   386.42, false),
    #                     "WB5"            => (2, 954.82,   1146.4, true),
    #                     "case6ww"        => (1, 2986,     2986,   false),
    #                     "case9"          => (2, 373.8,    1458.8, true),
    #                     "case9mod"       => (2, 234.6,    1320.4, true),
    #                     "case14"         => (2, 721.5,    5371.5, true),
    #                     # "case22loop"     => (1, 4538.8,   4538.8, false), ## Absent from data repo...
    #                     "case30"         => (2, 268.915,  316.49, true)]) #,
    #                     # "case39"         => (1, 1887.2,   1887.2, false)) #,
    #                     # "case39mod1"     => (2, 773.36,   942.34, true),
    #                     # "case39mod2"     => (1, 940.34,   940.34, false),
    #                     # "case57"         => (1, 25338,    25338,  false),
    #                     # "case89pegase"   => (1, 5817.6,   5817.6, false)) #,
    #                     # "case118"        => (1, 86298,    86298,  false),
    #                     # "case118mod"     => (1, 86079,    86079,  false),
    #                     # "case300"        => (2, 254840.4, 475430, true),
    #                     # "case300mod"     => (1, 287950,   287950, false),
    #                     # "case1354pegase" => (2, 74053,    74053,  true))

    testfolder = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "tests", "sos_example6")

    @testset "instance $instance, (rmeqs=$rmeqs)" for (instance, order) in Iterators.filter(x->x[2]==1, keys(OPFsols)), rmeqs in Set([false])

        problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
        problem = pb_cplx2real(problem_c)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = order,
                                            params = Dict(:opt_outlev=>0))

        primobj, dualobj = run_hierarchy(problem, relax_ctx);

        ε_rel = OPFsols[(instance, 1)].rel_opt_gap
        ε_abs = max(primobj, dualobj) * ε_rel

        @show primobj, dualobj,  OPFsols[(instance, 1)].primal_solvalue, ε_abs, instance

        @test primobj ≈ OPFsols[(instance, 1)].primal_solvalue atol=ε_abs
        @test dualobj ≈ primobj atol=ε_abs
    end
end
