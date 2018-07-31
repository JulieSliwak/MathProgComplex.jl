using ArgParse, MathProgComplex, OPFInstances

function parse_commandline()
  s = ArgParseSettings()
  s.description = "Load the *.sdp files at the given location and start Mosel solve."
  @add_arg_table s begin
    "instance_path"
        help=".dat instance path"
        arg_type = String
        required = true
    "d"
        help="relaxation order"
        arg_type = Int
        required = true
    "logpath"
        help = "folder storing all output"
        arg_type = String
        required = true
    "relaxation type"
        help = "either 0 (SOSRelaxation) or 1 (MomentRelaxation)"
        arg_type = Int
        required = true
  end
  return parse_args(s)
end


function main(args)
    input_params = parse_commandline()

    # input_params = Dict("instance_path"=>getinstancepath("Matpower", "QCQP", "case9"),
    #                     "d"=>1,
    #                     "logpath"=>joinpath("tltu", "tita"),
    #                     "relaxation type"=>0)

    instance_path = input_params["instance_path"]
    d = input_params["d"]
    hierarchykind = :Real
    symmetries = DataType[]
    logpath = input_params["logpath"]
    pbsolved = input_params["relaxation type"]==0 ? (:MomentRelaxation) : (:SOSRelaxation)

    instance_name = splitdir(instance_path)[2][1:end-4]
    mkpath(logpath)

    @printf("->  Launching %30s  %2d  %s\n", instance_name, d, logpath)

    ## Run functions once for precompilation...
    WB2_path = getinstancepath("Matpower", "QCQP", "WB2")
    problem_C, point = import_from_dat(WB2_path)
    problem = pb_cplx2real(problem_C)
    relax_ctx = set_relaxation(problem; hierarchykind=hierarchykind,
                                        d=1,
                                        params = Dict(:opt_outlev=>0,
                                                      :opt_outmode=>0,
                                                      :opt_solver=>:None))
    run_hierarchy(problem, relax_ctx; save_pbs=false)


    ## Build real problem
    problem_C, point = import_from_dat(instance_path)
    problem = pb_cplx2real(problem_C)

    ## Build relaxation context
    relax_ctx = set_relaxation(problem; hierarchykind=hierarchykind,
                                        d=d,
                                        symmetries=symmetries,
                                        params = Dict(:opt_outlev=>1,
                                                      :opt_outmode=>1,
                                                      :opt_outcsv=>1,
                                                    #   :opt_msk_maxtime=>2*3600,
                                                      :opt_pbsolved=>pbsolved,
                                                      :opt_solver=>:None,
                                                      :opt_outname=>joinpath(logpath, "momentsos.log"),
                                                      :opt_outcsvname=>joinpath(logpath, "momentsos_solve.csv")))
    relax_ctx.relaxparams[:pb_name] = instance_name

    run_hierarchy(problem, relax_ctx; save_pbs=false)

    final_output(relax_ctx)

end

main(ARGS)