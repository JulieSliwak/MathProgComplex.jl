using OPFInstances, BenchmarkTools, DataStructures

println("Loading module MathProgComplex")
tic();
using MathProgComplex
toc();

function main(instance = case30pwl)
    instancepath = getinstancepath("Matpower", "QCQP", instance)
    # instancepath = getinstancepath("Matpower", "QCQP", "case89pegase")
    # instancepath = getinstancepath("Matpower", "QCQP", "case300")
    println("\n\nWorking on $(splitdir(instancepath))")

    seek_efficiency!(false)
    @show seek_efficiency()
    @show MathProgComplex.DictType

    pb_c, pt = import_from_dat(instancepath)
    problem = pb_cplx2real(pb_c)

    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                        # symmetries=[PhaseInvariance],
                                        d = 1,
                                        params = Dict(:opt_outlev=>0,
                                                      :opt_outmode=>0,
                                                      :opt_outcsv=>0))

    max_cliques = get_maxcliques(relax_ctx, problem)
    relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)

    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

    println("--- time build_sparsity")
    @btime build_sparsity($relax_ctx, $problem, $max_cliques)

    mmtrel_pb = build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    println("--- time build_momentrelaxation")
    @btime build_momentrelaxation($relax_ctx, $problem, $momentmat_param, $localizingmat_param, $max_cliques)

    sdpinstance = build_SOSrelaxation(relax_ctx, mmtrel_pb)

    println("--- time build_SOSrelaxation")
    @btime build_SOSrelaxation($relax_ctx, $mmtrel_pb)

    path = joinpath(pwd(), "Mosek_runs", "worksdp")
    mkpath(path)
    export_SDP(sdpinstance, path)
    sdp_instance = read_SDPPrimal(path)

    println("--- time read_SDPPrimal")
    @btime read_SDPPrimal($path)


    # sdp = SDP_Problem()

    # set_constraints!(sdp, sdp_instance)
    # set_vartypes!(sdp, sdp_instance)
    # set_blocks!(sdp, sdp_instance)
    # set_linvars!(sdp, sdp_instance)

    # set_matrices!(sdp, sdp_instance)
    # set_linear!(sdp, sdp_instance)
    # set_const!(sdp, sdp_instance)

    # # println(sdp)

    # primal = SortedDict{Tuple{String,String,String}, Float64}()
    # dual = SortedDict{Tuple{String, String, String}, Float64}()

    # primobj, dualobj = solve_mosek(sdp::SDP_Problem, primal, dual, sol_info=relax_ctx.relaxparams)

    # # final_output(relax_ctx)

    return
end

function main_work()
    instancepath = getinstancepath("Matpower", "QCQP", "case300")
    println("\nWorking on $(splitdir(instancepath))")

    seek_efficiency!(true)
    @show seek_efficiency()
    @show MathProgComplex.DictType

    pb_c, pt = import_from_dat(instancepath)

    pb_cplx2real(pb_c)
    seek_efficiency!(false)

    return
end


# main_work()
main("WB2")
main("case30pwl")
main("case89pegase")
main("case300")
