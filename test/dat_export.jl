using MathProgComplex

@testset "dat export import consistency check" begin
    # min x - y
    # st  x + x^2 + x*y + y^2 <= 1
    #     -2 <= x, y <= 2
    # solution: x+y = -1/3
    # optimal objective -1-4/sqrt(3)

    qcqp = Problem()

    x = MathProgComplex.Variable("x", Real)
    y = MathProgComplex.Variable("y", Real)

    set_objective!(qcqp, x - y)

    add_constraint!(qcqp, "ctr", (x + x^2 + x*y + y^2) << 1)
    add_constraint!(qcqp, "x bounds", -2 << x << 2)
    add_constraint!(qcqp, "y bounds", -2 << y << 2)

    point = Point()
    dat_exportpath = joinpath(Pkg.dir("MathProgComplex"), "tmp_datexport", "myinstance.dat")
    mkpath(splitdir(dat_exportpath)[1])

    export_to_dat(qcqp, dat_exportpath, point)
    qcqp2 = import_from_dat(dat_exportpath)


    amplscriptpath = joinpath(Pkg.dir("MathProgComplex"), "src", "export_dat")
    # run_knitro(dat_exportpath, amplscriptpath)

    @test qcqp.objective == qcqp2.objective

    @test Set(keys(qcqp.constraints)) == Set(keys(qcqp2.constraints))
    for (ctrname, ctr) in qcqp.constraints
        @test qcqp.constraints[ctrname].ub == qcqp2.constraints[ctrname].ub
        @test qcqp.constraints[ctrname].lb == qcqp2.constraints[ctrname].lb
        @test qcqp.constraints[ctrname].p == qcqp2.constraints[ctrname].p
    end
end
