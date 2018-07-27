using Mosek

using MathProgComplex, DataStructures, OPFInstances

!isdefined(:MPC) && (const MPC = MathProgComplex)

# include(joinpath(Pkg.dir("MathProgComplex"), "src", "SDPhierarchy", "SDP_Instance", "build_from_SDPDual.jl"))
# include(joinpath(Pkg.dir("MathProgComplex"), "src", "SDPhierarchy", "io", "export_SDP_Instance.jl"))

function main()
    # problem_c, point = import_from_dat(getinstancepath("Matpower", "QCQP", "case9"))
    # problem = pb_cplx2real(problem_c)


    # order = 1
    problem_fct = lasserre_ex1

    problem, order_to_obj, ε_abs = problem_fct()

    testsolver= :MosekSolver
    pbsolved = :MomentRelaxation

    d = 2
    relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                d = d,
                                params = Dict(:opt_outlev=>0,
                                            :opt_solver=>testsolver,
                                            :opt_pbsolved=>pbsolved))

    primobj, dualobj = run_hierarchy(problem, relax_ctx)

    @show pbsolved, problem_fct, d
    @show d, primobj, dualobj, order_to_obj[d]

    @show primobj - order_to_obj[d]
    @show primobj - dualobj

    # @assert isapprox(primobj, order_to_obj[d], atol=1e-2)
    # @assert isapprox(dualobj, primobj, atol=1e-2)

    order = 2

    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                         d = order,
    #                         params = Dict(:opt_outlev=>1,
    #                                     :opt_solver=>testsolver,
    #                                     :opt_pbsolved=>:SOSRelaxation))

    # primobj_sos, dualobj_sos = run_hierarchy(problem, relax_ctx)

    # relax_ctx = set_relaxation(problem; hierarchykind=:Real,
    #                     d = order,
    #                     params = Dict(:opt_outlev=>1,
    #                                 :opt_solver=>testsolver,
    #                                 :opt_pbsolved=>:MomentRelaxation))

    # primobj_mmt, dualobj_mmt = run_hierarchy(problem, relax_ctx)

    # @show primobj_sos, dualobj_sos
    # @show primobj_mmt, dualobj_mmt
    # return

    workpath = "SDPwork"
    ispath(workpath) && rm(workpath, recursive = true)
    mkpath(workpath)

    worksos = joinpath(workpath, "sos"); mkpath(worksos)
    workmmt = joinpath(workpath, "mmt"); mkpath(workmmt)

    relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                            issparse = false,
                                            d = order,
                                            params = Dict(:opt_outlev=>1))


    max_cliques = MPC.get_maxcliques(relax_ctx, problem)
    relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)
    momentmat_param, localizingmat_param = MPC.build_sparsity(relax_ctx, problem, max_cliques)

    ########################################
    # Build the moment, SOS relaxation problem
    momentrel = MPC.build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)
    sosrel = MPC.build_SOSrelaxation(relax_ctx, momentrel)


    sdp_instance_moment::SDP_Problem = build_SDP_Instance_from_SDPDual(momentrel)
    sdp_instance_sos::SDP_Problem = MPC.build_SDP_Instance_from_SDPPrimal(sosrel)


    prim_sos = SortedDict{Tuple{String,String,String}, Float64}()
    dual_sos = SortedDict{Tuple{String, String, String}, Float64}()

    primobj_sos, dualobj_sos = MPC.solve_JuMP(sdp_instance_sos::MPC.SDP_Problem, :MosekSolver, prim_sos,
                                                                      dual_sos,
                                                                      sol_info=relax_ctx.relaxparams,
                                                                      logname = joinpath(worksos, "Mosek.log"),
                                                                      printlog = true,
                                                                      debug = false,
                                                                      optsense=:Max)

    # println(STDOUT, "----- SOS relaxation solutions -----")
    # println(STDOUT, "Dual solution")
    # for ((blockname, var1, var2), val) in dual_sos
    #     @printf(STDOUT, "%30s %10s %10s %f\n", blockname, var1, var2, val)
    # end

    # println(STDOUT, "Primal solution")
    # for ((blockname, var1, var2), val) in prim_sos
    #     @printf(STDOUT, "%30s %10s %10s %f\n", blockname, var1, var2, val)
    # end

    prim_mmt = SortedDict{Tuple{String, String, String}, Float64}()
    dual_mmt = SortedDict{Tuple{String, String, String}, Float64}()



    primobj_mmt, dualobj_mmt = MPC.solve_mosek(sdp_instance_moment::MPC.SDP_Problem, prim_mmt,
                                                                      dual_mmt,
                                                                      sol_info=relax_ctx.relaxparams,
                                                                      logname = joinpath(".", "Mosek.log"),
                                                                      printlog = true,
                                                                      debug = true,
                                                                      optsense=:Min)

    primobj_mmt, dualobj_mmt = MPC.solve_JuMP(sdp_instance_moment::MPC.SDP_Problem, :MosekSolver, prim_mmt,
                                                                    dual_mmt,
                                                                    sol_info=relax_ctx.relaxparams,
                                                                    logname = joinpath(workmmt, "Mosek.log"),
                                                                    printlog = true,
                                                                    debug = true,
                                                                    optsense = :Min)
    # println(STDOUT, "----- Moment relaxation solutions -----")
    # println(STDOUT, "Primal solution")
    # for ((blockname, var1, var2), val) in prim_mmt
    #     @printf(STDOUT, "%30s %10s %10s %f\n", blockname, var1, var2, val)
    # end

    # println(STDOUT, "Dual solution")
    # for ((blockname, var1, var2), val) in dual_mmt
    #     @printf(STDOUT, "%30s %10s %10s %f\n", blockname, var1, var2, val)
    # end

    # for k in keys(dual_mmt)
    #     dual_mmt[k] *= -1
    # end

    for k in keys(dual_sos)
        dual_sos[k] *= -1
    end

    # println(STDOUT, "----- Moment primal | SOS dual -----")
    # printcompared(prim_mmt, dual_sos)

    # println(STDOUT, "----- Moment dual | SOS primal -----")
    # printcompared(dual_mmt, prim_sos)

    # for (mmt_blockname, mmt_var1, mmt_var2) in keys(prim_mmt)
    #     @printf("%30s  %10s  %10s  %i\n", mmt_blockname, mmt_var1, mmt_var2, mmt_var1 >= mmt_var2)
    # end

    # for (mmt_blockname, mmt_var1, mmt_var2) in keys(prim_sos)
    #     @printf("%30s  %10s  %10s  %i\n", mmt_blockname, mmt_var1, mmt_var2, mmt_var1 >= mmt_var2)
    # end

    @printf("sos relaxation ; primal: %10f   dual : %10f        sol: %10f\n", primobj_sos, dualobj_sos, order_to_obj[order])
    @printf("mmt relaxation ; primal: %10f   dual : %10f        sol: %10f\n", primobj_mmt, dualobj_mmt, order_to_obj[order])

    # @show SortedSet(collect(keys(prim_mmt)))
    # @show SortedSet(collect(keys(dual_sos)))

    return
end

function printcompared(prim_mmt, dual_sos)
    @printf("%30s %10s %10s %10s  |  %10s  %10s %10s %30s ||     %10s \n", "mat Z_i", "row id j", "col id k", "Z_i[j,k]", "Z_i*[j*,k*]", "row id j*", "col id k*", "mat Z_i*", "rel. difference (%)")
    for (k, v) in zip(prim_mmt, dual_sos)
        (mmt_blockname, mmt_var1, mmt_var2), mmt_val = k
        (sos_blockname, sos_var1, sos_var2), sos_val = v

        reldiff = 100 * abs(mmt_val-sos_val) / max(abs(mmt_val), abs(sos_val), 1e-3)

        @printf("%30s %10s %10s % 10f  |  % 10f  %10s %10s %30s || delta is %10.2f \n", mmt_blockname, mmt_var1, mmt_var2, mmt_val, sos_val, sos_var1, sos_var2, sos_blockname, reldiff)
    end
end

main()