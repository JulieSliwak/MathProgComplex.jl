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
## Build polynomial problem
x1 = Variable("x", Real)
x2 = Variable("y", Real)
problem = Problem()
set_objective!(problem, -1.0*x1)
add_constraint!(problem, "ineq", (x1^2+x2^2) << 4)
θ1 = π/3
add_constraint!(problem, "eq_rot1", (cos(θ1)*x1+sin(θ1)*x2) == 0)

## Set parameters
relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                    d = 1,
                                    params = Dict(:opt_outlev=>1,
                                                  :opt_pbsolved=>:SOSRelaxation,
                                                  :opt_solver=>:MosekCAPI))

## Build the order 1 SOS relaxation and pass it to Mosek via the C API.
primobj, dualobj = run_hierarchy(problem, relax_ctx);
```

Set important options, get objective values, get solution values, set order, symmetries, sparse

## simple example

Here is the low level equivalent of the previous code example.

```julia
## Build polynomial problem
x1 = Variable("x", Real)
x2 = Variable("y", Real)
problem = Problem()
set_objective!(problem, -1.0*x1)
add_constraint!(problem, "ineq", (x1^2+x2^2) << 4)
θ1 = π/3
add_constraint!(problem, "eq_rot1", (cos(θ1)*x1+sin(θ1)*x2) == 0)

relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                        d = 1,
                                        params = Dict(:opt_outlev=>1))


########################################
# Build sparsity pattern, chordal extension, maximal cliques.
# Simple dense problem : one clique with all variables.
max_cliques = MPC.get_maxcliques(relax_ctx, problem)
relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)

########################################
# Compute moment and localizing matrices parameters: order et variables
momentmat_param, localizingmat_param = MPC.build_sparsity(relax_ctx, problem, max_cliques)

########################################
# Build the moment relaxation problem
momentrel = MPC.build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)
# Convert to a primal SDP problem

########################################
sosrel = MPC.build_SOSrelaxation(relax_ctx, momentrel)

# println("\n--------------------------------------------------------")
# println("sosrel = \n$sosrel")
# println("--------------------------------------------------------")

mkdir(joinpath(workpath, "SOSpb"))
MPC.export_SDPPrimal(sosrel, joinpath(workpath, "SOSpb"), renamemoments=false)

# sdp_instance = MPC.build_SDP_Instance_from_sdpfiles(path)

# sdp_instance_moment::SDP_Problem = build_SDP_Instance_from_SDPDual(momentrel)
sdp_instance_sos::MPC.SDP_Problem = MPC.build_SDP_Instance_from_SDPPrimal(sosrel)

# mkpath(joinpath(workpath, "export_sdp_pb_sos"))
# mkpath(joinpath(workpath, "export_sdp_pb_moment"))

# export_SDP_Instance(sdp_instance_sos, joinpath(workpath, "export_sdp_pb_sos"))
# export_SDP_Instance(sdp_instance_moment, joinpath(workpath, "export_sdp_pb_moment"))
# MPC.export_SDPPrimal(sosrel, joinpath(workpath, "export_sdp_pb_moment"), renamemoments=false)


# export_SDP_Instance(sdp_instance_moment, workpath)


primal = SortedDict{Tuple{String,String,String}, Float64}()
dual = SortedDict{Tuple{String, String, String}, Float64}()

primobj, dualobj = MPC.solve_JuMP(sdp_instance_sos, :CSDPSolver, primal, dual;
                                                            logname = "Mosek_run.log",
                                                            printlog = false,
                                                            msk_maxtime = relax_ctx.relaxparams[:opt_msk_maxtime],
                                                            sol_info = relax_ctx.relaxparams,
                                                            optsense = :Max)

```