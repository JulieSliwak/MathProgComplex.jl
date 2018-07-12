using MathProgComplex, DataStructures
include(joinpath(Pkg.dir("MathProgComplex"), "src", "SDPhierarchy", "SDP_Instance", "build_from_SDPDual.jl"))


function main()

    # problem = buildPOP_WB2(setnetworkphase=false)
    problem = buildPOP_WB2()

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = 1,
                                        params = Dict(:opt_outlev=>1,
                                                      :opt_outmode=>0,
                                                      :opt_outcsv=>0))

    problem, relax_ctx = lasserre_ex5(d=1)

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
    momentrel = build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    println("\n--------------------------------------------------------")
    println("momentrel = $momentrel")
    println("--------------------------------------------------------")

    ########################################
    # Convert to a primal SDP problem
    sosrel = build_SOSrelaxation(relax_ctx, momentrel)

    # println("\n--------------------------------------------------------")
    # println("sosrel = \n$sosrel")
    # println("--------------------------------------------------------")

    path = joinpath(pwd(), "Mosek_runs", "worksdp")
    mkpath(path)
    export_SDP(sosrel, path, renamemoments=false)


    sdp_instance = build_SDP_Instance_from_sdpfiles(path)

    # println("\n--------------------------------------------------------")
    # println(sdp_instance)
    # println("--------------------------------------------------------")

    sdp_instance = build_SDP_Instance_from_SDPPrimal(sosrel)

    # println("\n--------------------------------------------------------")
    # println(sdp_instance)
    # println("--------------------------------------------------------")

    sdp_instance = build_SDP_Instance_from_SDPDual(momentrel)

    println("\n--------------------------------------------------------")
    println(sdp_instance)
    println("--------------------------------------------------------")


    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp_instance::SDP_Problem, primal, dual, sol_info=relax_ctx.relaxparams,
                                                                            debug=true)

    # final_output(relax_ctx)
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
