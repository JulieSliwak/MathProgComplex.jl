export import_from_dat

"""
    pb, point = import_from_dat(instancepath::String, precondcstrspath::String)

Build the polynomial optimization problem described by the `instancepath` file,
with the dispensory preconditionning descritpion `precondcstrspath` along with
the initial point possibly provided in the file (defaults value is null).
"""
function import_from_dat(instancepath::String; filename::String="real_minlp_instance.dat", precondfilename::String="")
    point = Point()
    variables = SortedDict{String, Variable}()
    exponents = SortedDict{String, Exponent}()
    pb = Problem()

    instance_str = open(joinpath(instancepath, filename))
    l = jump_comments!(instance_str)


    ## Collect and define variables
    line = matchall(r"\S+", l)
    while line[1] == "VAR_TYPE" && !eof(instance_str)
        if line[2] == "REAL"
            var = Variable(line[3], Real)
        elseif line[2] == "BOOL"
            var = Variable(line[3], Bool)
        elseif line[2] == "C"
            var = Variable(line[3], Complex)
        else
            error("import_to_dat(): Unknown variable type $(line[2]) for variable $(line[3]).")
        end

        variables[line[3]] = var
        add_variable!(pb, var)

        setindex!(point, convert(var.kind, parse(line[5]) + im*parse(line[6])), var)

        ## Mark where objective definition begins
        mark(instance_str)

        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end


    ## Move forward to the monomial definition section
    while line[1] != "MONO_DEF" && !eof(instance_str)
        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end

    ## Build all monomials
    while line[1] == "MONO_DEF" && !eof(instance_str)
        exponame = line[2]
        var = variables[line[3]]
        if !haskey(exponents, exponame)
            exponents[exponame] = Exponent()
        end
        exponents[exponame] = product(exponents[exponame], Exponent(SortedDict(var=>Degree(parse(line[5]), parse(line[6])))))
        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end

    ## Add final monomial
    if line[1] == "MONO_DEF"
        exponame = line[2]
        var = variables[line[3]]
        if !haskey(exponents, exponame)
            exponents[exponame] = Exponent()
        end
        exponents[exponame] = product(exponents[exponame], Exponent(SortedDict(var=>Degree(parse(line[5]), parse(line[6])))))
    end

    ## Reset stream to the objective defintion
    reset(instance_str)
    l = readline(instance_str)
    line = matchall(r"\S+", l)

    ## Build polynomial objective
    p = Polynomial()
    while line[2] == "OBJ"
        λ = parse_λ(line[5], line[6])
        var1, var2 = line[3:4]
        if line[1] == "MONO"
            p += λ * exponents[var1]
        else
            p += λ * (var1!="NONE" ? conj(variables[var1]) : 1) * (var2!="NONE" ? variables[var2] : 1)
        end
        l = readline(instance_str)
        line = matchall(r"\S+", l)
    end
    set_objective!(pb, p)

    ## Build constraints
    next_state = :AssembleCtr
    cur_ctr = line[2]
    var1, var2 = line[3:4]

    lb = -Inf
    ub = +Inf
    p = Polynomial()
    λ = parse_λ(line[5], line[6])
    while next_state!=:ReadLine || (!eof(instance_str) && line[1] != "MONO_DEF")

        state = next_state
        if state == :ReadLine
            ## Readline
            l = readline(instance_str)
            line = matchall(r"\S+", l)

            var1, var2 = line[3:4]
            λ = parse_λ(line[5], line[6])

            ## Determine whether current ctr should be completed or new ctr should be set
            if line[2] == cur_ctr
                next_state = :AssembleCtr
            else
                next_state = :SaveCtr
            end

            ## If monomial difnition section reached, exit loop
            if line[1] == "MONO_DEF"
                next_state = :ReadLine 
            end

        elseif state == :AssembleCtr
            if line[1] == "MONO"
                p += λ * exponents[var1]
            elseif line[1] ∈ SortedSet(["CONST", "LIN", "QUAD"])
                p += λ * (var1!="NONE" ? conj(variables[var1]) : 1) * (var2!="NONE" ? variables[var2] : 1)
            elseif line[1] == "UB"
                ub = λ
            elseif line[1] == "LB"
                lb = λ
            else
                error("import_from_dat(): Unknown variable type $(line[1]) for constraint $(line[2]).")
            end

            next_state = :ReadLine

        elseif state == :SaveCtr
            add_constraint!(pb, String(cur_ctr), lb << p << ub)

            lb = -Inf
            ub = +Inf
            p = Polynomial()
            cur_ctr = line[2]
            next_state = :AssembleCtr
        else
            error("import_from_dat(): Unknown state $state")
        end
    end

    # Add final constraint
    add_constraint!(pb, String(cur_ctr), lb << p << ub)

    ## Set preconditioning flag
    if precondfilename != ""
        precond_str = open(joinpath(instancepath, precondfilename))

        l = jump_comments!(precond_str)
        line = matchall(r"\S+", l)
        while !eof(precond_str)
            if line[2] == "SQRT"
                pb.constraints[line[1]].precond = :sqrt
            else
                error("import_from_dat(): Unknown preconditioning $(line[2]) for constraint $(line[1]).")
            end
            l = readline(precond_str)
            line = matchall(r"\S+", l)
        end
    end

    return pb, point
end


"""
    l = jump_comments!(io::IOStream)

Jump comments, i.e. lines starting by '#', in the `io` stream, and return the
first non commented line.
"""
function jump_comments!(io::IOStream)
    l = ""
    while !eof(io)
        l = readline(io)
        ismatch(r"\s*#", l) || break
    end
    return l
end

function parse_λ(realpart::T, imagpart::T) where T <: Union{String, SubString}
    λ = 0
    if realpart == "Inf"
        λ += Inf
    elseif realpart == "-Inf"
        λ += -Inf
    else
        λ += parse(realpart)
    end

    if imagpart == "Inf"
        λ += Inf*im
    elseif imagpart == "-Inf"
        λ += -Inf*im
    else
        λ += parse(imagpart)*im
    end
    return λ
end 