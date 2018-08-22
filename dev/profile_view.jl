using MathProgComplex, ProfileView

function main()

    instance_path = getinstancepath("Matpower", "QCQP", "case9")
    ## Build real problem
    problem_C, point = import_from_dat(instance_path)
    problem = pb_cplx2real(problem_C)

    pbsolved = :MomentRelaxation

    ## Build relaxation context
    relax_ctx = set_relaxation(problem; hierarchykind=hierarchykind,
                                        d=d,
                                        symmetries=symmetries,
                                        params = Dict(:opt_outlev=>1,
                                                      :opt_outmode=>0,
                                                      :opt_outcsv=>0,
                                                    #   :opt_solver_maxtime=>2*3600,
                                                      :opt_relaxationkind=>pbsolved,
                                                      :opt_solver=>:None,
                                                      :opt_outname=>joinpath(logpath, "momentsos.log"),
                                                      :opt_outcsvname=>joinpath(logpath, "momentsos_solve.csv")))

    max_cliques = get_maxcliques(relax_ctx, problem)

    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

    mmtrel_pb = build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques);

    sdp_pb = build_SDP_Instance_from_SDPDual(mmtrel_pb);

    Profile.clear()
    @profile build_momentrelaxation(relax_ctx, problem, momentmat_param, localizingmat_param, max_cliques)
    ProfileView.view()

end