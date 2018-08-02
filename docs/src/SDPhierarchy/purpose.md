# Guide

[TODO]

## Purpose - general principle of relaxation

- SDP problem
- relaxing POPs...
- to SDP problems of increasing size, converging to POP objective

- Case of OPFs, POP arising in managing AC high level current transportation network: order 2 is fine up to hundreds of variables, order 1 is not always enough... Hence the will to find outs in the middle.

## Features of the hierarchy

- real or complex,
  - real has been known since 2001, applied to real POPs
  - complex is usefull when input POP are complex. Good to avoid converting CPOP to RPOP as this op. looses structural information
  - difficulty: complex SDP solvers don't exist yet
- multiordered,
  - a relaxation order is defined for each constraint
  - those will define the size of the problem
- dense or sparse,
  - decompose one large problem into several smaller problems and coupling constraints
  - valid decomposition is not easy to derive, depends on the constraints and relaxation orders
- symmetries
  - symmetries present in objective and all constraints of the POP are respected by non-null moments

## references

## simple "high level" workflow

```julia
x1 = Variable("x", Real)
x2 = Variable("y", Real)
problem = Problem()
set_objective!(problem, -1.0*x1)
add_constraint!(problem, "ineq", (x1^2+x2^2) << 4)
θ1 = π/3
add_constraint!(problem, "eq_rot1", (cos(θ1)*x1+sin(θ1)*x2) == 0)

relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                    d = 1,
                                    params = Dict(:opt_outlev=>1,
                                                  :opt_pbsolved=>:MomentRelaxation))

primobj, dualobj = run_hierarchy(problem, relax_ctx);
```

Set important options, get objective values, get solution values, set order, symmetries, sparse

## simple example

In depth example here: ....