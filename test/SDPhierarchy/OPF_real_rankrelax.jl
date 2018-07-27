using MathProgComplex, OPFInstances
using Base.Test


brokentestset = Set([])

testcases = OrderedSet([("WB2", 1),
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

@testset "SDP hierarchy ACOPFs relaxation - rank 1" begin

    @testset "$pbsolved tests" for pbsolved in [:SOSRelaxation, :MomentRelaxation]
        @testset "$instance problem - rank $d" for (instance, d) in testcases
            problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", instance))
            problem = pb_cplx2real(problem_c)

            relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                                d = d,
                                                params = Dict(:opt_outlev=>0,
                                                                :opt_solver=>testsolver,
                                                                :opt_pbsolved=>pbsolved))

            primobj, dualobj = run_hierarchy(problem, relax_ctx);

            ε_rel = OPFsols[(instance, 1)].rel_opt_gap
            ε_abs = max(primobj, dualobj) * ε_rel

            @show testsolver, pbsolved, instance, d
            @show primobj, dualobj,  OPFsols[(instance, 1)].primal_solvalue, ε_abs, instance

            @test primobj ≈ OPFsols[(instance, 1)].primal_solvalue atol=ε_abs
            @test dualobj ≈ primobj atol=ε_abs
        end
    end
end