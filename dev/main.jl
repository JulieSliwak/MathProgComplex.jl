using DataStructures, SCS, Mosek, OPFInstances

import MathProgComplex

!isdefined(:MPC) && (const MPC = MathProgComplex)

function main()
    # problem = buildPOP_WB2(setnetworkphase=false)
    # problem_c, pt = MPC.import_from_dat(getinstancepath("Matpower", "QCQP", "case9"))

    # problem = MPC.pb_cplx2real(problem_c)
    
    workpath = joinpath("Mosek_runs", "worksdp")
    ispath(workpath) && rm(workpath, recursive=true)
    mkpath(workpath)
    
    problem, relax = MPC.lasserre_ex1()
    
    cstobj = problem.objective[MPC.Exponent()]
    cstobj = 0
    problem.objective += -cstobj
    # @assert !haskey(problem.objective, MPC.Exponent())
    
    # problem_c = buildPOP_1v1c()
    # problem = MPC.Problem()
    # x1 = MPC.Variable("x", Real)
    # x2 = MPC.Variable("y", Real)
    
    # MPC.set_objective!(problem, -x1 + 1 - 0.14159)
    # MPC.add_constraint!(problem, "ctr2", (x1^4+x2^4) << 2^4)
    # MPC.add_constraint!(problem, "ctr1", (x1^2+x2^2) << 10^2)
    # println(problem)

    relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                            issparse = false,
                                            d = 2,
                                            params = Dict(:opt_outlev=>1,
                                                          :opt_outmode=>2,
                                                          :opt_outcsv=>0,
                                                          :opt_relaxationkind=>:SOSRelaxation,
                                                          :opt_solver=>:MosekSolver,
                                                          :opt_debug=>true))

    primobj_JuMP, dualobj_JuMP = MPC.run_hierarchy(problem, relax_ctx, indentedprint=true, save_pbs=true)


    relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                                issparse = false,
                                                d = 2,
                                                params = Dict(:opt_outlev=>1,
                                                            :opt_outmode=>2,
                                                            :opt_outcsv=>0,
                                                            :opt_relaxationkind=>:SOSRelaxation,
                                                            :opt_solver=>:MosekCAPI,
                                                            :opt_debug=>true))

    primobj_CAPI, dualobj_CAPI = MPC.run_hierarchy(problem, relax_ctx, indentedprint=true, save_pbs=true)

    @show primobj_JuMP, dualobj_JuMP

    @show primobj_CAPI, dualobj_CAPI

    @show cstobj
    return nothing

    # problem, relax_ctx = lasserre_ex1()
    # problem, relax_ctx = lasserre_ex2()
    # problem, relax_ctx = lasserre_ex3()
    # problem, relax_ctx = lasserre_ex5(d=2)
    # export_to_dat(problem, workpath, filename="QCQP.dat", exportprecond=false)

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    # max_cliques = MPC.get_case9cliques(relax_ctx, problem)
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
                                                                msk_maxtime = relax_ctx.relaxparams[:opt_solver_maxtime],
                                                                sol_info = relax_ctx.relaxparams,
                                                                optsense = :Max)

    @show primobj, dualobj

    return nothing
    # primobj, dualobj = MPC.solve_mosek(sdp_instance_sos::MPC.SDP_Problem, primal,
    #                                                                   dual,
    #                                                                   sol_info=relax_ctx.relaxparams,
    #                                                                   debug=false,
    #                                                                   logname = joinpath("Mosek.log"),
    #                                                                   optsense=:Max)


    # open(joinpath(workpath, "SOSrelsol.dat"), "w") do fout
    #     println(fout, "----- SOS relaxation solutions -----")
    #     println(fout, "Primal solution")
    #     for ((blockname, var1, var2), val) in primal
    #     @printf(fout, "%30s %10s %10s %f\n", blockname, var1, var2, val)
    #     end

    #     println(fout, "Dual solution")
    #     for ((blockname, var1, var2), val) in dual
    #     @printf(fout, "%30s %10s %10s %f\n", blockname, var1, var2, -val)
    #     end
    # end

    # primobjmom, dualobjmom = MPC.solve_mosek(sdp_instance_moment::MPC.SDP_Problem, primal,
    #                                                                   dual,
    #                                                                   sol_info=relax_ctx.relaxparams,
    #                                                                   debug=true,
    #                                                                   optsense=:Min)

    # open(joinpath(workpath, "mmtrelsol.dat"), "w") do fout
    #     println(fout, "----- Moment relaxation solutions -----")
    #     println(fout, "Dual solution")
    #     for ((blockname, var1, var2), val) in dual
    #     @printf(fout, "%30s %10s %10s %f\n", blockname, var1, var2, -val)
    #     end

    #     println(fout, "Primal solution")
    #     for ((blockname, var1, var2), val) in primal
    #     @printf(fout, "%30s %10s %10s %f\n", blockname, var1, var2, val)
    #     end
    # end

    # println("SOS relaxation:    $primobj, $dualobj")
    # println("Moment relaxation: $primobjmom, $dualobjmom")


    println("Objectives : $primobj, $dualobj")

    moseksolver = MosekSolver()

    # moseksolver = MosekSolver(LOG = 0)
    m = JuMP_from_SDP_Problem(sdp_instance_sos, moseksolver)

    JuMP.solve(m)

    return m

    objective = JuMP.getobjectivevalue(m)

    println()
    println("primobj    ", primobj)
    println("dualobj    ", dualobj)
    println("objective  ", objective)

    return problem
end

main()
