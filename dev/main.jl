using MathProgComplex, DataStructures

function main()

    # problem = buildPOP_WB2(v2max=1.022, setnetworkphase=false)
    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                                     # symmetries=[PhaseInvariance],
    #                                     d = 1)

    # problem = buildPOP_WB2(setnetworkphase=false)
    problem = buildPOP_WB5()
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = 2,
                                        params = Dict(:opt_outlev=>1,
                                                      :opt_outmode=>0,
                                                      :opt_outcsv=>0))

    # println("\n--------------------------------------------------------")
    # println("problem = \n$problem")

    # println("\n--------------------------------------------------------")
    # println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    max_cliques = get_maxcliques(relax_ctx, problem)
    relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)

    # println("\n--------------------------------------------------------")
    # println("max cliques =")
    # println(max_cliques)

    ########################################
    # Compute moment and localizing matrices parameters: order et variables
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

    # println("\n--------------------------------------------------------")
    # println("moment params =")
    # for (cliquename, dcl) in momentmat_param
    #     println("Moment matrix, $cliquename \t -> dcl = $dcl")
    # end
    # for (key, (val1, val2)) in localizingmat_param
    #     @printf("%15s \t -> di-ki = %i, \tcliques = ", key, val2)
    #     for clique in val1 print("$clique, ") end
    #     @printf("\b\b \n")
    # end

    ########################################
    # Build the moment relaxation problem
    mmtrel_pb = build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    # println("\n--------------------------------------------------------")
    # println("mmtrel_pb = $mmtrel_pb")

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)

    # println("\n--------------------------------------------------------")
    # println("sdpinstance = \n$sdpinstance")

    path = joinpath(pwd(), "Mosek_runs", "worksdp")
    mkpath(path)
    export_SDP(sdpinstance, path)
    sdp_instance = read_SDPPrimal(path)

    # println("VAR_TYPES size:     $(size(sdp_instance.VAR_TYPES))")
    # println("BLOCKS size:        $(size(sdp_instance.BLOCKS))")
    # println("LINEAR size:        $(size(sdp_instance.LINEAR))")
    # println("CONST size:         $(size(sdp_instance.CONST))")

    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_vartypes!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)
    set_linvars!(sdp, sdp_instance)

    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)

    # println(sdp)

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual, sol_info=relax_ctx.relaxparams)

    final_output(relax_ctx)
    # # println("Primal solution")
    # # for ((blockname, var1, var2), val) in primal
    # # @printf("%15s %5s %5s %f\n", blockname, var1, var2, val)
    # # end

    # # println("\nDual solution NEGATED")
    # # for var in problem.variables
    # #     ctrname = get_momentcstrname()
    # #     var1 = var[1]
    # #     var2 = "1"
    # #     val = dual[(ctrname, var1, var2)]
    # #     println("($(ctrname), $(var1), $(var2)) = $(-val)")
    # # end

    # println("Objectives : $primobj, $dualobj")
end

main()
