using Mosek

using MathProgComplex, DataStructures

!isdefined(:MPC) && (const MPC = MathProgComplex)

# include(joinpath(Pkg.dir("MathProgComplex"), "src", "SDPhierarchy", "SDP_Instance", "build_from_SDPDual.jl"))
# include(joinpath(Pkg.dir("MathProgComplex"), "src", "SDPhierarchy", "io", "export_SDP_Instance.jl"))

function main()
    problem = Problem()

    x1 = Variable("x1", Real); add_variable!(problem, x1)
    x2 = Variable("x2", Real); add_variable!(problem, x2)
    x3 = Variable("x3", Real); add_variable!(problem, x3)

    set_objective!(problem, -2*x1+x2-x3)
    add_constraint!(problem, "ctr1", (x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24) >> 0)
    add_constraint!(problem, "ctr2", (x1+x2+x3) << 4)
    add_constraint!(problem, "ctr3", (3*x2+x3) << 6)
    add_constraint!(problem, "def_x1", 0 << x1 << 2)
    add_constraint!(problem, "def_x2", 0 << x2)
    add_constraint!(problem, "def_x3", 0 << x3 << 3)


    problem = Problem()
    x1 = Variable("x", Real); add_variable!(problem, x1)
    x2 = Variable("y", Real); add_variable!(problem, x2)

    set_objective!(problem, -x1)
    # add_constraint!(problem, "ctr2", (x1^4+x2^4) << 2^4)
    add_constraint!(problem, "ctr1", (x1^2+x2^2) << 2^2)

    # problem = buildPOP_WB2(setnetworkphase=false)

    order_to_obj = SortedDict(1=>-6.0000,
                            2=>-5.6923,
                            3=>-4.0685,
                            4=>-4.0000)

    order = 1

    workpath = "SDPwork"
    ispath(workpath) && rm(workpath, recursive = true)
    mkpath(workpath)

    worksos = joinpath(workpath, "sos"); mkpath(worksos)
    workmmt = joinpath(workpath, "mmt"); mkpath(workmmt)

    relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                            issparse = false,
                                            d = 1,
                                            params = Dict(:opt_outlev=>1))


    max_cliques = MPC.get_maxcliques(relax_ctx, problem)
    relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)
    momentmat_param, localizingmat_param = MPC.build_sparsity(relax_ctx, problem, max_cliques)

    ########################################
    # Build the moment, SOS relaxation problem
    momentrel = MPC.build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)
    sosrel = MPC.build_SOSrelaxation(relax_ctx, momentrel)

    println(sosrel)


    sdp_instance_moment::SDP_Problem = build_SDP_Instance_from_SDPDual(momentrel)
    sdp_instance_sos::SDP_Problem = MPC.build_SDP_Instance_from_SDPPrimal(sosrel)


    prim_sos = SortedDict{Tuple{String,String,String}, Float64}()
    dual_sos = SortedDict{Tuple{String, String, String}, Float64}()

    primobj_sos, dualobj_sos = MPC.solve_mosek(sdp_instance_sos::MPC.SDP_Problem, prim_sos,
                                                                      dual_sos,
                                                                      sol_info=relax_ctx.relaxparams,
                                                                      logname = joinpath(worksos, "Mosek.log"),
                                                                      printlog = true,
                                                                      debug=true,
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
                                                                      logname = joinpath(workmmt, "Mosek.log"),
                                                                      printlog = true,
                                                                      debug = true,
                                                                      optsense=:Min)


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

    # for k in keys(dual_sos)
    #     dual_sos[k] *= -1
    # end

    println(STDOUT, "----- Moment primal | SOS dual -----")
    printcompared(prim_mmt, dual_sos)

    println(STDOUT, "----- Moment dual | SOS primal -----")
    printcompared(dual_mmt, prim_sos)

    # for (mmt_blockname, mmt_var1, mmt_var2) in keys(prim_mmt)
    #     @printf("%30s  %10s  %10s  %i\n", mmt_blockname, mmt_var1, mmt_var2, mmt_var1 >= mmt_var2)
    # end

    # for (mmt_blockname, mmt_var1, mmt_var2) in keys(prim_sos)
    #     @printf("%30s  %10s  %10s  %i\n", mmt_blockname, mmt_var1, mmt_var2, mmt_var1 >= mmt_var2)
    # end

    @printf("sos relaxation ; primal: %10f   dual : %10f        sol: %10f\n", primobj_sos, dualobj_sos, order_to_obj[order])
    @printf("mmt relaxation ; primal: %10f   dual : %10f        sol: %10f\n", primobj_mmt, dualobj_mmt, order_to_obj[order])

    @show SortedSet(collect(keys(prim_mmt)))
    @show SortedSet(collect(keys(dual_sos)))

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