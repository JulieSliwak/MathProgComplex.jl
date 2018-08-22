using MathProgComplex
using Base.Test

problems = [lasserre_ex1, lasserre_ex2, lasserre_ex5, test_gloptipoly_ex1, test_gloptipoly_ex2]

brokentestset = Set([(:MomentRelaxation, lasserre_ex2, 3), (:MomentRelaxation, lasserre_ex5, 2)])

@testset "SDP hierarchy base tests" begin

    @testset "$pbsolved tests" for pbsolved in [:SOSRelaxation, :MomentRelaxation]
        @testset "$problem_fct problem" for problem_fct in problems
            problem, order_to_obj, ε_abs = problem_fct()

            @testset "Order $d" for d in keys(order_to_obj)
                relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = d,
                                            params = Dict(:opt_outlev=>0,
                                                        :opt_solver=>testsolver,
                                                        :opt_relaxationkind=>pbsolved))

                primobj, dualobj = run_hierarchy(problem, relax_ctx)

                @show pbsolved, problem_fct, d
                @show d, primobj, dualobj, order_to_obj[d]

                if (pbsolved, problem_fct, d) in brokentestset
                    @test_broken primobj ≈ order_to_obj[d] atol=ε_abs
                else
                    @test primobj ≈ order_to_obj[d] atol=ε_abs
                    @test dualobj ≈ primobj atol=ε_abs*min(abs(primobj), abs(dualobj))
                end
            end
        end
    end
end