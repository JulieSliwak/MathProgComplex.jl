export solve_JuMP, JuMP_from_SDP_Problem

# for complex to real conversion, see klorel/mathprogcomplex a78593860b662fc1c09dbd3910656b84425729f7

"""
  (primalobj, dualobj) = solve_JuMP(problem, solver, primal, dual; debug, logname, printlog, msk_maxtime, sol_info, optsense)

Calls any JuMP interfaced SDP solver on `problem::SDP_Problem`. Returns the primal and dual objectives if possible.

TODO: update arguments
*Arguments* :
- `problem::SDP_Problem`
- `primal::SortedDict{Tuple{String,String,String}, Float64}`: primal solution `x`,
- `dual::SortedDict{Tuple{String, String, String}, Float64}` : dual solution `s`,
- `debug`: default is false, dump Mosek loaded problem
- `printlog` : if true (default), Mosek will write its log to the console,
- `logname` : if specified, Mosek log will be written to the given file name,
- `msk_maxtime` : Mosek max computation time, in seconds. Default -1 means no limit,
- `sol_info` : Information on problem and solution status upon termination of solve,
- `optsense` : default is `:Max`.

**Note**:
- Mosek expects lower triangular terms of the coefficient matrices. Hence diagonal or non-diagonal terms will not be scaled.
"""
function solve_JuMP(problem::SDP_Problem, solver::T,
                                            primal::SortedDict{Tuple{String,String,String}, Float64},
                                            dual::SortedDict{Tuple{String, String, String}, Float64};
                                            debug = false,
                                            logname = "",
                                            printlog = true,
                                            msk_maxtime = -1,            # Default -1 means no time limit
                                            sol_info = OrderedDict(),
                                            optsense = :Max) where T<:MathProgBase.AbstractMathProgSolver

    empty!(primal)
    empty!(dual)

    m = JuMP_from_SDP_Problem(problem, solver)

    JuMP.setobjectivesense(m, optsense)

    JuMP.solve(m)

    debug && MathProgBase.writeproblem(m.internalModel, "myfile.jtask")
    debug && mv("myfile.jtask", "solve_JuMP_task.json", remove_destination = true)

    primobj = JuMP.getobjectivevalue(m)

    ## TODO: - build up primal, dual solutions...
    ##       - have functions signature converge

    return primobj, primobj
end

function solve_JuMP(problem::SDP_Problem, solver::Symbol,
                                            primal::SortedDict{Tuple{String,String,String}, Float64},
                                            dual::SortedDict{Tuple{String, String, String}, Float64};
                                            debug = false,
                                            logname = "",
                                            printlog = true,
                                            msk_maxtime = -1,            # Default -1 means no time limit
                                            sol_info = OrderedDict(),
                                            optsense = :Max)

    empty!(primal)
    empty!(dual)

    @assert solver in OrderedSet([:MosekSolver, :SCSSolver, :CSDPSolver])
    if solver == :MosekSolver
        options = Any[]

        !printlog && push!(options, (:MSK_IPAR_LOG, 0))

        mysolver = Mosek.MosekSolver(options)

    elseif solver == :SCSSolver
        options = Any[]

        !printlog && push!(options, (:verbose, 0))

        mysolver = SCS.SCSSolver(options)

    elseif solver == :CSDPSolver
        options = Dict{Symbol, Any}()

        !printlog && (options[:printlevel] = 0)

        mysolver = CSDP.CSDPSolver(options)
    end

    m = JuMP_from_SDP_Problem(problem, mysolver)

    JuMP.setobjectivesense(m, optsense)

    JuMP.solve(m)

    debug && MathProgBase.writeproblem(m.internalModel, "myfile.jtask")
    debug && mv("myfile.jtask", "solve_JuMP_task.json", remove_destination = true)

    primobj = JuMP.getobjectivevalue(m)

    ## TODO: - build up primal, dual solutions...
    ##       - have functions signature converge

    return primobj, primobj
end



"""
    m = JuMP_from_SDP_Problem(sdp_pb::SDP_Problem, mysolver::T) where T<:MathProgBase.AbstractMathProgSolver

Given a JuMP interfaced SDP solver, convert the input `SDP_Problem` into a JuMP model `m`.
"""
function JuMP_from_SDP_Problem(sdp_pb::SDP_Problem, mysolver::T) where T<:MathProgBase.AbstractMathProgSolver
    m = JuMP.Model(solver = mysolver)

    ## Variables
    Zi = Dict{String, Array{JuMP.Variable, 2}}()
    for (blockname, block) in sdp_pb.name_to_sdpblock
        n = length(block.var_to_id)
        var = JuMP.@variable(m, [1:n,1:n], SDP, basename=blockname)
        Zi[blockname] = var
    end


    varlin = Dict()
    for varname in keys(sdp_pb.scalvar_to_id)
        varlin[varname] = JuMP.@variable(m, basename=varname)
    end


    ## Constraints

    # Gathering constraint and objective keys
    ctr_names = SortedSet([k[1] for k in keys(sdp_pb.matrices)])
    union!(ctr_names, [k[1] for k in keys(sdp_pb.linear)])
    union!(ctr_names, [k for k in keys(sdp_pb.cst_ctr)])
    ctr_names = SortedSet(collect(ctr_names))

    objective_keys = SortedSet()
    for ctr_name in ctr_names
        if ctr_name[1] == "1" && ctr_name[2] == "1"
            delete!(ctr_names, ctr_name)
            push!(objective_keys, ctr_name)
        end
    end

    obj = 0
    bodys = Dict{Tuple{String, String, String}, JuMP.GenericAffExpr{Float64,JuMP.Variable}}([ctr_name=>0 for ctr_name in ctr_names])

    # Building constraints and objective bodys
    for ((ctr_name, blockname, var1, var2), f_αβ) in sdp_pb.matrices
        block = sdp_pb.name_to_sdpblock[blockname]
        n = length(block.var_to_id)
        i = block.var_to_id[var1]
        j = block.var_to_id[var2]

        Ai = sparse([i], [j], [1], n, n)
        if (i != j) #&& (ctr_name ∉ objective_keys)
            Ai = sparse([i; j], [j; i], [1.; 1.], n, n)
        end

        if ctr_name in ctr_names
            bodys[ctr_name] += vecdot( Zi[blockname], Ai ) * f_αβ
        else
            @assert ctr_name in objective_keys
            obj += vecdot( Zi[blockname], Ai ) * f_αβ
        end
    end

    for ((ctr_name, varname), f_αβ) in sdp_pb.linear
        !haskey(bodys, ctr_name) && (bodys[ctr_name] = 0)

        if ctr_name in ctr_names
            bodys[ctr_name] += varlin[varname] * f_αβ
        else
            @assert ctr_name in objective_keys
            obj += varlin[varname] * f_αβ
        end
    end

    for (ctr_name, f_αβ) in sdp_pb.cst_ctr
        if ctr_name in ctr_names
            bodys[ctr_name] += f_αβ
        else
            @assert ctr_name in objective_keys
            obj += f_αβ
        end
    end

    # Adding constraints and objective to JuMP model
    JuMP.@objective(m, Max, obj)
    JuMP.@constraintref ctrs[1:length(ctr_names)]

    for (i, ctr_name) in enumerate(ctr_names)
        ctrs[i] = JuMP.@constraint(m, bodys[ctr_name] == 0)
    end

    return m

end
