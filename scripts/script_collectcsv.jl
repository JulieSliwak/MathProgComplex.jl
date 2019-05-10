using CSV, DataFrames, ArgParse

# function parse_commandline()
#   s = ArgParseSettings()
#   s.description = "Merge all atomic csv into one large csv."
#   @add_arg_table s begin
#     "run_dir"
#     help = "Directory which contains each run folder."
#     arg_type = String
#     required = true
#   end
#   return parse_args(s)
# end

function main(args)
    println("--- Merge_csv.jl ---")

    ## Setting working directory
    # input_params = parse_commandline()
    # input_params = Dict("run_dir" => joinpath("Mosek_runs", "pararuns", "07_31-19h20"))
    input_params = Dict("run_dir" => joinpath("07_31-19h20"))
    workdir = input_params["run_dir"]

    info("Working in $workdir")

    run_folders = Set(filter(x->isdir(joinpath(workdir, x)), readdir(workdir)))
    println("Found $(length(run_folders)) folders")

    global_csv = DataFrames.DataFrame()

    n_addedinstances = 1
    for folder in run_folders
        curfolder = joinpath(workdir, folder)
        println("-- working on $curfolder")


        ## Collecting information
        csvfile = joinpath(curfolder, "momentsos_solve.csv")


        out_str = filter(x->ismatch(r"hierarchy.*\.o", x), readdir(curfolder))
        @assert length(out_str) == 1
        out_str = first(out_str)

        err_str = filter(x->ismatch(r"hierarchy.*\.e", x), readdir(curfolder))
        @assert length(err_str) == 1
        err_str = first(err_str)

        instancename = out_str[11:end-2]

        waskilled = readstring(joinpath(curfolder, out_str)) == "Killed\n"
        waserror = readstring(joinpath(curfolder, err_str)) != ""

        # # Has timed out ?
        # hastimedout = true
        # if isfile(csvfile) != 0
        #     hastimedout = false
        # end
        @assert isfile(csvfile)

        if isfile(csvfile)
            ## Reading csv
            cur_csv = CSV.read(csvfile, delim=";")

            cur_csv[:slv_run_status] = "All went well."
            if waskilled
                cur_csv[:slv_run_status] = "Process was killed (by qsub)."
                info("Process was killed (qsub)")
            end
            if waserror
                cur_csv[:slv_run_status] = "Process errored (likely was qdel)."
                info("Process errored (likely was qdel).")
            end

            if eltype(cur_csv[:pb_name]) == Missing
                warn("Empty pb name")
                newcsv = DataFrame()
                for name in names(global_csv)
                    if name == :pb_name
                        newcsv[:pb_name] = instancename
                    else
                        newcsv[name] = cur_csv[name]
                    end
                end
                cur_csv = newcsv
            end

            allowmissing!(cur_csv)

            if size(global_csv) == (0,0)
                global_csv = deepcopy(cur_csv)
                allowmissing!(global_csv)
            elseif names(global_csv) != names(cur_csv)
                # @show names(cur_csv)
                # @show length(names(cur_csv))
                # @show names(global_csv)
                # @show length(names(global_csv))
                warn("Ignoring  $folder :   Names")
            elseif DataFrames.eltypes(global_csv) != DataFrames.eltypes(cur_csv)
                # @show DataFrames.eltypes(global_csv)
                # @show DataFrames.eltypes(cur_csv)
                # @show names cur_csv
                warn("Ignoring  $folder :   Elttypes")
                # append!(global_csv, cur_csv)
            else
                n_addedinstances += 1
                append!(global_csv, cur_csv)
            end
        end

        println("   Size global CSV : $(size(global_csv))")
    end

    csvname = "report_log.csv"
    println("Writing full csv file $(joinpath(workdir, csvname))")

    CSV.write(joinpath(workdir, csvname), global_csv, delim=';')

    nbfine = size(global_csv[global_csv[:slv_run_status].=="All went well.", :], 1)

    info("$nbfine runs were fine.")
    info("$n_addedinstances / $(length(run_folders)) added to final csv.")

    return nothing
end

main(ARGS)
