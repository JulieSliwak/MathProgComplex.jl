include("..\\src\\MathProgComplex.jl")
using Test

using Ipopt

rosenbrock = Problem()

x = Variable("x", Real)
y = Variable("y", Real)
set_objective!(rosenbrock, Base.power_by_squaring((1-x),2) + 100*Base.power_by_squaring((y-Base.power_by_squaring(x,2)),2))

mysolver = IpoptSolver(print_level = 0)
m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(rosenbrock, mysolver)
print(m)
solve(m)

sol = get_JuMP_solution(m, variables_jump, rosenbrock)


@test sol[x] ≈ 1.0 atol=1e-5
@test sol[y] ≈ 1.0 atol=1e-5
@test get_objective(rosenbrock, sol) ≈ 0.0 atol=1e-5
