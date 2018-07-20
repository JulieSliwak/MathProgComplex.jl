using DataStructures, SCS, Mosek, OPFInstances

import MathProgComplex
import JuMP

!isdefined(:MPC) && (const MPC = MathProgComplex)

# include(joinpath(Pkg.dir("MathProgComplex"), "src", "SDPhierarchy", "SDP_Instance", "build_from_SDPDual.jl"))
include(joinpath(Pkg.dir("MathProgComplex"), "src", "SDPhierarchy", "solvers", "JuMP.jl"))

function main()
    # problem = buildPOP_WB2(setnetworkphase=false)

    workpath = joinpath("Mosek_runs", "worksdp")
    ispath(workpath) && rm(workpath, recursive=true)
    mkpath(workpath)

    problem_c, pt = MPC.import_from_dat(getinstancepath("Matpower", "QCQP", "WB5"))
    problem = MPC.pb_cplx2real(problem_c)

    println(problem)

    relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                        issparse = false,
                                        d = 1,
                                        params = Dict(:opt_outlev=>1,
                                                      :opt_outmode=>0,
                                                      :opt_outcsv=>0))

    # problem, relax_ctx = MPC.lasserre_ex5(d=1)

    # println("\n--------------------------------------------------------")
    # println("problem = \n$problem")

    # println("\n--------------------------------------------------------")
    # println("relax_ctx = \n$relax_ctx")

    ########################################
    # Construction du sparsity pattern, extension chordale, cliques maximales.
    # max_cliques = MPC.get_case9cliques(relax_ctx, problem)
    max_cliques = MPC.get_maxcliques(relax_ctx, problem)

    relax_ctx.relaxparams[:opt_nb_cliques] = length(max_cliques)

    ########################################
    # Compute moment and localizing matrices parameters: order et variables
    momentmat_param, localizingmat_param = MPC.build_sparsity(relax_ctx, problem, max_cliques)

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
    momentrel = MPC.build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)

    # println("\n--------------------------------------------------------")
    # println("momentrel = $momentrel")
    # println("--------------------------------------------------------")

    ########################################
    # Convert to a primal SDP problem
    sosrel = MPC.build_SOSrelaxation(relax_ctx, momentrel)

    # println("\n--------------------------------------------------------")
    # println("sosrel = \n$sosrel")
    # println("--------------------------------------------------------")

    # path = joinpath(pwd(), "Mosek_runs", "worksdp")
    # mkpath(path)

    # sdp_instance = MPC.build_SDP_Instance_from_sdpfiles(path)

    sdp_instance = MPC.build_SDP_Instance_from_SDPPrimal(sosrel)
    MPC.export_SDPPrimal(sosrel, workpath, renamemoments=false)

    # sdp_instance = build_SDP_Instance_from_SDPDual(momentrel)

    # println("\n--------------------------------------------------------")
    # println(sdp_instance)
    # println("--------------------------------------------------------")


    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj, dualobj = MPC.solve_mosek(sdp_instance::MPC.SDP_Problem, primal,
                                                                      dual,
                                                                      sol_info=relax_ctx.relaxparams,
                                                                      debug=false)

    # final_output(relax_ctx)

    # println("Primal solution")
    # for ((blockname, var1, var2), val) in primal
    # @printf("%15s %10s %10s %f\n", blockname, var1, var2, val)
    # end

    # println("Dual solution")
    # for ((blockname, var1, var2), val) in dual
    # @printf("%15s %10s %10s %f\n", blockname, var1, var2, -val)
    # end

    println("Objectives : $primobj, $dualobj")

    moseksolver = MosekSolver()
    m = JuMP_from_SDP_Problem(sdp_instance, moseksolver)

    JuMP.solve(m)

    objective = JuMP.getobjectivevalue(m)

    println()
    println("primobj    ", primobj)
    println("dualobj    ", dualobj)
    println("objective  ", objective)
    return 
end

main()
