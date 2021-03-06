using DataStructures, OPFInstances

function main()

    #############################################################
    # Parameters
    nmax_order1 = 350
    nmax_order2 = 82
    juliapath = joinpath("~", "software", "julia-0.6.3", "bin", "julia")
    juliascript = joinpath(pwd(), "dev", "script_dat_to_hierarchysimple.jl")
    qsubwkdir = pwd()

    date = String(Dates.format(now(), "mm_dd-HHhMM"))
    workdir = joinpath(pwd(), "Mosek_runs", "pararuns", date)

    ispath(workdir) && rm(workdir, recursive=true)
    mkpath(workdir)


    params = OrderedSet()
    tempfolders = Set()

    matpower_path = first(splitdir(getinstancepath("Matpower", "QCQP", "WB2")))
    instances = OrderedSet(sort(readdir(matpower_path), by=x->parse(first(matchall(r"\d+", x)))))

    instances_order1 = filter(x->parse(first(matchall(r"\d+", x))) ≤ nmax_order1, instances)
    instances_order2 = filter(x->parse(first(matchall(r"\d+", x))) ≤ nmax_order2, instances)


    for instance in instances_order1
        tempfolder = randstring(4)
        while tempfolder in tempfolders
            tempfolder = randstring(4)
        end
        push!(tempfolders, tempfolder)

        instancename = first(split(instance, "."))
        push!(params, (instancename, joinpath(matpower_path, instance), 1, joinpath(workdir, tempfolder), 0))
    end

    for instance in instances_order1
        tempfolder = randstring(4)
        while tempfolder in tempfolders
            tempfolder = randstring(4)
        end
        push!(tempfolders, tempfolder)

        instancename = first(split(instance, "."))
        push!(params, (instancename, joinpath(matpower_path, instance), 1, joinpath(workdir, tempfolder), 1))
    end

    for instance in instances_order2
        tempfolder = randstring(4)
        while tempfolder in tempfolders
            tempfolder = randstring(4)
        end
        push!(tempfolders, tempfolder)

        instancename = first(split(instance, "."))
        push!(params, (instancename, joinpath(matpower_path, instance), 2, joinpath(workdir, tempfolder), 0))
    end

    for instance in instances_order2
        tempfolder = randstring(4)
        while tempfolder in tempfolders
            tempfolder = randstring(4)
        end
        push!(tempfolders, tempfolder)

        instancename = first(split(instance, "."))
        push!(params, (instancename, joinpath(matpower_path, instance), 2, joinpath(workdir, tempfolder), 1))
    end

    for param in params
        println(param)
    end

    sleep(3)

    for param in params
        instancename, instancepath, d, logpath, hierarchykind = param
        # println(`qsub -e $logpath/hierarchy_$(instancename).e -o $logpath/hierarchy_$(instancename).o $juliapath dev/script_dat_to_hierarchysimple.jl $instancename $d $logpath`)
        # run(`qsub -e $logpath/hierarchy_$(instancename).e -o $logpath/hierarchy_$(instancename).o -b y julia dev/script_dat_to_hierarchysimple.jl $instancepath $d $logpath`)
        mkpath(logpath)
        open(joinpath(logpath, "runscript.sh"), "w") do f
            println(f, "cd $(pwd())")
            println(f, "qsub -e $logpath/hierarchy_$(instancename).e -o $logpath/hierarchy_$(instancename).o -wd $qsubwkdir -b y $juliapath $juliascript $instancepath $d $logpath $hierarchykind")
        end

        run(`bash $(joinpath(logpath, "runscript.sh"))`)
    end

end

main()