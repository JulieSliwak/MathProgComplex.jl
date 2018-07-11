using Base.Test
using MathProgComplex

function run_pb_orders(problem, order_to_obj)
    data = SortedDict()
    for d in keys(order_to_obj)
        relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                            d = d,
                                            params = Dict(:opt_outlev=>0))

        logpath = joinpath(Pkg.dir("MathProgComplex"), "Mosek_runs", "Lasserreex", "d_$d")
        ispath(logpath) && rm(logpath, recursive=true)
        mkpath(logpath)
        cur_obj, dualobj = run_hierarchy(problem, relax_ctx, logpath, indentedprint=true, save_pbs=true)

        @show d, cur_obj, dualobj
        data[d] = (cur_obj, dualobj)
    end
    return data
end

@testset "GloptiPoly small real problems" begin

    ## http://homepages.laas.fr/henrion/papers/gpcocos.pdf p7
    @testset "example1" begin
        problem = Problem()
        x1 = Variable("x1", Real); add_variable!(problem, x1)
        x2 = Variable("x2", Real); add_variable!(problem, x2)
        x3 = Variable("x3", Real); add_variable!(problem, x3)

        set_objective!(problem, -2*x1+x2-x3)
        add_constraint!(problem, "ctr1", (x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24) >> 0)
        add_constraint!(problem, "ctr2", (x1+x2+x3) << 4)
        add_constraint!(problem, "ctr3", (3*x2+x3) << 6)
        add_constraint!(problem, "def_x1", 0 << x1 << 2)
        add_constraint!(problem, "def_x2", 0 << x2)
        add_constraint!(problem, "def_x3", 0 << x3 << 3)

        ε=1e-3
        order_to_obj = SortedDict(1=>-6.0000,
                                2=>-5.6923,
                                3=>-4.0685,
                                4=>-4.0000)

        data_sol = run_pb_orders(problem, order_to_obj)

        @printf("%5s  %15s  %15s  | %15s", "order", "primal obj", "dual obj", "expected obj\n")
        for (d, (primobj, dualobj)) in data_sol
            @printf("%5i  %15f  %15f  | %15f\n", d, primobj, dualobj, order_to_obj[d])
            @test primobj ≈ order_to_obj[d] atol=ε
            @test dualobj ≈ primobj atol=mosek_optgap*min(abs(primobj), abs(dualobj))
        end
    end

    @testset "example2" begin
            ## http://homepages.laas.fr/henrion/papers/gloptipoly.pdf p14
        problem = Problem()
        x1 = Variable("x1", Real); add_variable!(problem, x1)
        x2 = Variable("x2", Real); add_variable!(problem, x2)
        x3 = Variable("x3", Real); add_variable!(problem, x3)
        x4 = Variable("x4", Real); add_variable!(problem, x4)
        x5 = Variable("x5", Real); add_variable!(problem, x5)

        set_objective!(problem, 42*x1+44*x2+45*x3+47*x4+47.5*x5-50*(x1^2+x2^2+x3^2+x4^2+x5^2))
        add_constraint!(problem, "ctr1", (20*x1+12*x2+11*x3+7*x4+4*x5) << 40)
        add_constraint!(problem, "def_x1", 0 << x1 << 1)
        add_constraint!(problem, "def_x2", 0 << x2 << 1)
        add_constraint!(problem, "def_x3", 0 << x3 << 1)
        add_constraint!(problem, "def_x4", 0 << x4 << 1)
        add_constraint!(problem, "def_x5", 0 << x5 << 1)

        ε=1e-3
        order_to_obj = SortedDict(1=>NaN,         # order 1 relaxation is primal infeasable
                                2=>-17.9189,
                                3=>-17.0000)


        data_sol = run_pb_orders(problem, order_to_obj)

        @printf("%5s  %15s  %15s  | %15s", "order", "primal obj", "dual obj", "expected obj\n")
        for (d, (primobj, dualobj)) in data_sol
            @printf("%5i  %15f  %15f  | %15f\n", d, primobj, dualobj, order_to_obj[d])
            if d==1
                @test isnan(primobj)
                @test isnan(dualobj)
            else
                @test primobj ≈ order_to_obj[d] atol=ε
                @test dualobj ≈ primobj atol = mosek_optgap*min(abs(primobj), abs(dualobj))
            end
        end
    end
end