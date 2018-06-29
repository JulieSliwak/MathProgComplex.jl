
function run_hierarchy(problem::Problem, relax_ctx::RelaxationContext, logpath; indentedprint=false,
                                                                                max_cliques::Dict{String, Set{Variable}}=Dict{String, Set{Variable}}(),
                                                                                save_pbs=false)

    open(joinpath(logpath, "pb_opt.log"), "w") do fout
        print(fout, problem)
    end

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    if max_cliques == Dict{String, Set{Variable}}()
        max_cliques = get_maxcliques(relax_ctx, problem)
    end
    relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)

    ########################################
    # Compute moment matrices parameters: order et variables
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

    ########################################
    # Compute moment and localization matrices
    mmtrel_pb, t, bytes, gctime, memallocs = @timed MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques);
    # mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)
    relax_ctx.relaxparams[:slv_mmtrel_t] = t
    relax_ctx.relaxparams[:slv_mmtrel_bytes] = bytes

    if save_pbs
        open(joinpath(logpath, "mmt_pb.log"), "w") do fout
            print(fout, mmtrel_pb)
        end
    end

    ########################################
    # Convert to a primal SDP problem
    sdpinstance, t, bytes, gctime, memallocs = @timed build_SOSrelaxation(relax_ctx, mmtrel_pb);
    # sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)
    relax_ctx.relaxparams[:slv_sosrel_t] = t
    relax_ctx.relaxparams[:slv_sosrel_bytes] = bytes

    export_SDP(sdpinstance, logpath, indentedprint=indentedprint)

    sdp, t, bytes, gctime, memallocs = @timed build_mosekpb(logpath);
    relax_ctx.relaxparams[:slv_mskstruct_t] = t
    relax_ctx.relaxparams[:slv_mskstruct_bytes] = bytes

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj = dualobj = NaN
    try
        if (relax_ctx.relaxparams[:opt_outmode]!=1) && (relax_ctx.relaxparams[:opt_outlev] ≥ 1)
            printlog = true
        else
            printlog = false
        end

        primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual; logname = joinpath(logpath, "Mosek_run.log"),
                                                                       printlog = printlog,
                                                                       msk_maxtime = relax_ctx.relaxparams[:opt_msk_maxtime],
                                                                       sol_info = relax_ctx.relaxparams)

    catch err
        relax_ctx.relaxparams[:slv_prosta] = err.msg
    end

    params_file = joinpath(logpath, "maxcliques_relaxctx.txt")
    isfile(params_file) && rm(params_file)
    open(params_file, "w") do fcliques
        println(fcliques, "# max_cliques are:")
        println(fcliques, max_cliques)
        println(fcliques, "# relaxation_ctx is:")
        println(fcliques, relax_ctx)
    end

    return primobj, dualobj
end

function build_relaxation(problem::Problem, relax_ctx::RelaxationContext; max_cliques::Dict{String, Set{Variable}} = Dict{String, Set{Variable}}())

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    if max_cliques == Dict{String, Set{Variable}}()
        max_cliques = get_maxcliques(relax_ctx, problem)
    end

    ########################################
    # Compute moment matrices parameters: order et variables
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

    ########################################
    # Compute moment and localization matrices
    mmtrel_pb = MomentRelaxation{Float64}(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    ########################################
    # Convert to a primal SDP problem
    sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)

    return sdpinstance

end


"""
    sdp = build_mosekpb(logpath)

    Build the sdp problem described at `logpath` into the generic structure interfaced with SDP solvers (Mosek for now).
"""
function build_mosekpb(logpath::String)
    sdp_instance = read_SDPInstance(logpath)

    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_vartypes!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)
    set_linvars!(sdp, sdp_instance)

    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)
    return sdp
end


function build_mosekpb(SOS_pb::SDPInstance, logpath::String; indentedprint=false)
    export_SDP(SOS_pb, logpath, indentedprint=indentedprint)

    sdp = build_mosekpb(logpath)

    sdp_instance = read_SDPInstance(logpath)

    sdp = SDP_Problem()

    set_constraints!(sdp, sdp_instance)
    set_vartypes!(sdp, sdp_instance)
    set_blocks!(sdp, sdp_instance)
    set_linvars!(sdp, sdp_instance)

    set_matrices!(sdp, sdp_instance)
    set_linear!(sdp, sdp_instance)
    set_const!(sdp, sdp_instance)
    return sdp
end

