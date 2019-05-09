using JuMP, Ipopt

@testset "dat export import consistency check" begin
    # min x - y
    # st  x + x^2 + x*y + y^2 <= 1
    #     -2 <= x, y <= 2
    # solution: x+y = -1/3
    # optimal objective -1-4/sqrt(3)

    qcqp = Problem()

    x = MPC.Variable("x", Complex)
    y = MPC.Variable("y", Real)
    b = MPC.Variable("b", Bool)

    set_objective!(qcqp, x - y)

    add_constraint!(qcqp, "ctr1", (x + x^2 + x*y + y^2) << 1)
    add_constraint!(qcqp, "ctr2", (x^3 + (1.0+3im)*conj(x^2)*x*y + y^2) << (1 + 2im))
    add_constraint!(qcqp, "ctr3", -3im << (5*b + 5im*y) << 1)
    add_constraint!(qcqp, "x_bounds", -2 << x << 2)
    add_constraint!(qcqp, "y_bounds", -2 << y << 2)

    point = Point([x], [3])
    dat_exportpath = joinpath(pwd(),"tmp_datexport")
    # dat_exportpath = joinpath(Pkg.dir("MathProgComplex"), "tmp_datexport")
    mkpath(splitdir(dat_exportpath)[1])

    export_to_dat(qcqp, dat_exportpath, filename="POP.dat", point=point)

    qcqp2, point2 = import_from_dat(joinpath(dat_exportpath, "POP.dat"))

    rm(joinpath(dat_exportpath, "POP.dat"))

    @test qcqp.variables == qcqp2.variables
    @test qcqp.objective == qcqp2.objective

    @test Set(keys(qcqp.constraints)) == Set(keys(qcqp2.constraints))

    for (ctrname, ctr) in qcqp.constraints
        @test qcqp.constraints[ctrname].ub == qcqp2.constraints[ctrname].ub
        @test qcqp.constraints[ctrname].lb == qcqp2.constraints[ctrname].lb
        @test qcqp.constraints[ctrname].p == qcqp2.constraints[ctrname].p
    end

    @test point == point2

    rm(dat_exportpath, recursive=true)
end

@testset "WB5.dat import and Ipopt solve" begin
    instancepath = joinpath(pwd(), "test", "instances")
    # instancepath = joinpath(Pkg.dir("MathProgComplex"), "test", "instances")
    WB5_cplx, initpt = import_from_dat(joinpath(instancepath, "WB5.dat"))

    WB5 = pb_cplx2real(WB5_cplx)

    mysolver = IpoptSolver(print_level = 0)
    m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(WB5, mysolver)
    solve(m)

    sol = get_JuMP_solution(m, variables_jump, WB5)

    @test get_objective(WB5, sol) ≈ 1146.478 atol=1e-3
end
