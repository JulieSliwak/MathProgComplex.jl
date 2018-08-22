using DataStructures, SCS, Mosek, OPFInstances
using MathProgComplex

using JuMP, KNITRO

!isdefined(:MPC) && (const MPC = MathProgComplex)

function main()
    problem_c, pt = MPC.import_from_dat(getinstancepath("Matpower", "QCQP", "WB5"))
    problem = MPC.pb_cplx2real(problem_c)

    workpath = joinpath("Mosek_runs", "worksdp")
    ispath(workpath) && rm(workpath, recursive=true)
    mkpath(workpath)

    # cstobj = problem.objective[MPC.Exponent()]
    # cstobj = 0
    # problem.objective += -cstobj
    # @assert !haskey(problem.objective, MPC.Exponent())

    di = Dict("LOAD_2_Im" => 1,
            "LOAD_2_Re" => 2,
            "LOAD_3_Im" => 1,
            "LOAD_3_Re" => 2,
            "LOAD_4_Im" => 1,
            "LOAD_4_Re" => 1,
            "UNIT_1_Im" => 1,
            "UNIT_1_Re" => 1,
            "UNIT_5_Im" => 1,
            "UNIT_5_Re" => 1,
            "VOLTM_1_Re" => 1,
            "VOLTM_2_Re" => 1,
            "VOLTM_3_Re" => 1,
            "VOLTM_4_Re" => 1,
            "VOLTM_5_Re" => 1)

    relax_ctx = MPC.set_relaxation(problem; hierarchykind=:Real,
                                            issparse = false,
                                            # d = 2,
                                            di = di,
                                            params = Dict(:opt_outlev=>1,
                                                          :opt_outmode=>0,
                                                          :opt_relaxationkind=>:MomentRelaxation,
                                                          :opt_solver=>:MosekCAPI))

    primal = SortedDict{Tuple{String,String,String}, Float64}()
    dual = SortedDict{Tuple{String, String, String}, Float64}()

    primobj_JuMP, dualobj_JuMP = MPC.run_hierarchy(problem, relax_ctx, indentedprint=true, save_pbs=true,
                                                                        primsol = primal,
                                                                        dualsol = dual)

    varname2ind = Dict([varname=>ind for (ind, varname) in enumerate(keys(problem.variables))])
    ind2varname = Dict([ind=>varname for (varname, ind) in varname2ind])


    n = length(varname2ind)

    mat = Matrix{Float64}(n, n)

    for ((clique, var1, var2), val) in primal
        if haskey(varname2ind, var1) && haskey(varname2ind, var2)
            mat[varname2ind[var1], varname2ind[var2]] = val
            mat[varname2ind[var2], varname2ind[var1]] = val
        end
    end

    eigenvals, eigenvecs = eig(mat)

    # mat2 = zeros(n, n)
    # for i=1:n
    #     mat2 += eigenvals[i] * (eigenvecs[: , i] * transpose(eigenvecs[: , i]))
    # end

    # println("\nEigenvals:")
    # display(eigenvals)

    # println("\nWorking with:")
    maxvec = eigenvecs[:, n]
    display(maxvec)

    ## Deal with sign
    if count(x->x<=0, maxvec) > n/2
        maxvec *= -1
    end

    ## Find the closest corresponding measure
    point_maxvec = Point()
    for i=1:n
        varname = ind2varname[i]
        vartype = problem.variables[varname]
        setindex!(point_maxvec, maxvec[i], Variable(varname, vartype))
    end


    localpb = Problem()
    localpb.variables = problem.variables

    obj = Polynomial()
    for ((clique, var1, var2), val) in primal
        if clique == "clique1" && var2 <= var1
            add!(obj, (parsevar(var1) * parsevar(var2) - val)^2)
        end
    end

    set_objective!(localpb, obj)

    solver = KnitroSolver()
    m, JuMPvar = get_JuMP_cartesian_model(localpb, solver)
    solve(m)

    pointknitro = Point()
    for (varname, jumpvar) in JuMPvar
        println("$jumpvar       ", getvalue(jumpvar))
        setindex!(pointknitro, getvalue(jumpvar), Variable(varname, Real))
    end

    println("--------------------------------------")
    println("\nslacks points - maxvec:\n", get_slacks(problem, point_maxvec))
    println("\nslacks points - knitro:\n", get_slacks(problem, pointknitro))

    return nothing
end

function parsevar(var::String)
    if var == "1"
        return 1.0
    elseif ismatch(r"\AVOLT_\d_(Re|Im)\Z", var)
        return Variable(var, Real)
    elseif ismatch(r"\AVOLT_\d_(Re|Im)\*VOLT_\d_(Re|Im)\Z", var)
        var1, var2 = matchall(r"VOLT_\d_(Re|Im)", var)
        return Variable(var1, Real) * Variable(var2, Real)
    elseif ismatch(r"\AVOLT_\d_(Re|Im)\^2\Z", var)
        var1 = matchall(r"VOLT_\d_(Re|Im)", var)
        return Variable(string(var1), Real)^2
    else
        @show var
        error()
    end
end

main()
