using JuMP, Ipopt

@testset "Rosenbrock Ipopt" begin
    # min (1-x)^2 + 100*(y-x^2)^2
    # solution: (x,y) = (1,1)
    # optimal objective 0
    rosenbrock = Problem()

    x = MPC.Variable("x", Real)
    y = MPC.Variable("y", Real)
    set_objective!(rosenbrock, (1-x)^2 + 100*(y-x^2)^2)

    mysolver = IpoptSolver(print_level = 0)
    m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(rosenbrock, mysolver)
    solve(m)

    sol = get_JuMP_solution(m, variables_jump, rosenbrock)

    @test sol[x] ≈ 1.0 atol=1e-5
    @test sol[y] ≈ 1.0 atol=1e-5
    @test get_objective(rosenbrock, sol) ≈ 0.0 atol=1e-5
end


@testset "QCQP Ipopt" begin
    # min x - y
    # st  x + x^2 + x*y + y^2 <= 1
    #     -2 <= x, y <= 2
    # solution: x+y = -1/3
    # optimal objective -1-4/sqrt(3)

    qcqp = Problem()

    x = MPC.Variable("x", Real)
    y = MPC.Variable("y", Real)

    set_objective!(qcqp, x - y)

    add_constraint!(qcqp, "ctr", (x + x^2 + x*y + y^2) << 1)
    add_constraint!(qcqp, "x bounds", -2 << x << 2)
    add_constraint!(qcqp, "y bounds", -2 << y << 2)

    mysolver = IpoptSolver(print_level = 0)
    m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(qcqp, mysolver)
    solve(m)

    sol = get_JuMP_solution(m, variables_jump, qcqp)

    ## Test for optimal solution
    @test sol[x] + sol[y] ≈ -1/3 atol=1e-5
    @test get_objective(qcqp, sol) ≈ -1-4/sqrt(3) atol=1e-5

    ## Test slack functions
    slacks = get_slacks(qcqp, sol)
    @test real(slacks.coords[MPC.Variable("x bounds", Complex)]) > 0
    @test real(slacks.coords[MPC.Variable("y bounds", Complex)]) > 0

    @test imag(slacks.coords[MPC.Variable("x bounds", Complex)]) == 0
    @test imag(slacks.coords[MPC.Variable("y bounds", Complex)]) == 0

    @test get_minslack(qcqp, sol)[1] < 1e-6
end


@testset "Matpower case9 Ipopt" begin
    instancepath = joinpath(Pkg.dir("MathProgComplex"), "test", "instances")
    pb_real, ~ = import_from_dat(joinpath(instancepath,"case9real.dat"))

    mysolver = IpoptSolver(print_level = 0)
    m, variables_jump = get_JuMP_cartesian_model(pb_real, mysolver)
    solve(m)
    @test getobjectivevalue(m) ≈ 1.45883471040144e+003 atol=1e-6
end
