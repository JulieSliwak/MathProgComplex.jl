using Base.Test, MathProgComplex, OPFInstances

## Beware, rank relaxation value is for problem without objective constant term

@testset "ACOPFs - real formulation - order 1" begin

    testsCases = OrderedSet([("WB2", 1),
                            ("LMBM3", 1),
                            ("WB5", 1),
                            ("case6ww", 1),
                            ("case9mod", 1),
                            ("case9", 1),
                            ("case14", 1)])
                            # ("case24_ieee_rts", 1),
                            # ("case_ieee30", 1),
                            # ("case30pwl", 1),
                            # ("case30", 1)])

    testfolder = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "tests", "sos_example6")

    @testset "instance $instance, (rmeqs=$rmeqs)" for (instance, order) in testsCases, rmeqs in Set([false])

        problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
        problem = pb_cplx2real(problem_c)

        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = order,
                                            params = Dict(:opt_outlev=>0,
                                                          :opt_solver=>testsSolver))

        primobj, dualobj = run_hierarchy(problem, relax_ctx);

        ε_rel = OPFsols[(instance, 1)].rel_opt_gap
        ε_abs = max(primobj, dualobj) * ε_rel

        @show primobj, dualobj,  OPFsols[(instance, 1)].primal_solvalue, ε_abs, instance

        @test primobj ≈ OPFsols[(instance, 1)].primal_solvalue atol=ε_abs
        @test dualobj ≈ primobj atol=ε_abs
    end
end
