using MathProgComplex, DataStructures, OPFInstances

function main()

    problem_c, pt = import_from_dat(getinstancepath("Matpower", "QCQP", "case9"))
    problem = pb_cplx2real(problem_c)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = 1,
                                        params = Dict(:opt_outlev=>1))

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
    sosrel = build_SOSrelaxation(relax_ctx, momentrel)

    # println("\n--------------------------------------------------------")
    # println("sosrel = \n$sosrel")


    path = joinpath(pwd(), "Mosek_runs", "worksdp")
    mkpath(path)
    export_SDP(sosrel, path)

    sdp_instance = build_SDP_Instance_from_sdpfiles(path)


    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = solve_mosek(sdp_instance::SDP_Problem, primal, dual, sol_info=relax_ctx.relaxparams)

end

main()
