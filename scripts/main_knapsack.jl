using MathProgComplex, DataStructures, OPFInstances

function main()

    problem = buildPOP_knapsack()

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        d = 1,
                                        binvar_d = 1,
                                        params = Dict(:opt_outlev=>1,
                                                      :opt_pbsolved=>:MomentRelaxation,
                                                      :opt_solver=>:MosekCAPI))


    # return run_hierarchy(problem, relax_ctx)
    # println("\n--------------------------------------------------------")
    # println("problem = \n$problem")

    # println("\n--------------------------------------------------------")
    # println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    max_cliques = get_maxcliques(relax_ctx, problem)
    relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)

    # println("\n--------------------------------------------------------")
    # println("max cliques =\nmax_cliques")

    ########################################
    # Compute moment and localizing matrices parameters: order et variables
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)


    ########################################
    # Build the moment relaxation problem
    momentrel = build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    # println("\n--------------------------------------------------------")
    # println("momentrel = $momentrel")

    ########################################
    # Convert to a primal SDP problem
    # sosrel = build_SOSrelaxation(relax_ctx, momentrel)

    # println("\n--------------------------------------------------------")
    # println("sosrel = \n$sosrel")

    sdp_instance = build_SDP_Instance_from_SDPDual(momentrel)
    # sdp_instance = build_SDP_Instance_from_SDPPrimal(sosrel)


    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp_instance::SDP_Problem, primal, dual, sol_info=relax_ctx.relaxparams, optsense=:Min)

    for i=1:5
        println("x$i:   ", primal[("clique1", "x$i", "1")])
    end

    for i=1:5
        println("x$i^2:   ", primal[("clique1", "x$i", "x$i")])
    end


    return primal
end

main()
