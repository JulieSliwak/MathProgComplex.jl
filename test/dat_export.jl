using MathProgComplex

@testset "dat export import consistency check" begin
    # min x - y
    # st  x + x^2 + x*y + y^2 <= 1
    #     -2 <= x, y <= 2
    # solution: x+y = -1/3
    # optimal objective -1-4/sqrt(3)

    qcqp = Problem()

    x = MathProgComplex.Variable("x", Complex)
    y = MathProgComplex.Variable("y", Real)
    b = MathProgComplex.Variable("b", Bool)

    set_objective!(qcqp, x - y)

    add_constraint!(qcqp, "ctr1", (x + x^2 + x*y + y^2) << 1)
    add_constraint!(qcqp, "ctr2", (x^3 + (1.0+3im)*conj(x^2)*x*y + y^2) << (1 + 2im))
    add_constraint!(qcqp, "ctr3", -3im << (5*b + 5im*y) << 1)
    add_constraint!(qcqp, "x_bounds", -2 << x << 2)
    add_constraint!(qcqp, "y_bounds", -2 << y << 2)

    point = Point([x], [3])
    dat_exportpath = joinpath(Pkg.dir("MathProgComplex"), "tmp_datexport")
    mkpath(splitdir(dat_exportpath)[1])

    export_to_dat(qcqp, dat_exportpath, filename="POP.dat", point=point)
    
    qcqp2, point2 = import_from_dat(dat_exportpath, filename="POP.dat")


    # amplscriptpath = joinpath(Pkg.dir("MathProgComplex"), "src", "export_dat")
    # run_knitro(dat_exportpath, amplscriptpath)

    @test qcqp.variables == qcqp2.variables
    @test qcqp.objective == qcqp2.objective

    @test Set(keys(qcqp.constraints)) == Set(keys(qcqp2.constraints))

    for (ctrname, ctr) in qcqp.constraints
        @test qcqp.constraints[ctrname].ub == qcqp2.constraints[ctrname].ub
        @test qcqp.constraints[ctrname].lb == qcqp2.constraints[ctrname].lb
        @test qcqp.constraints[ctrname].p == qcqp2.constraints[ctrname].p
    end
end
