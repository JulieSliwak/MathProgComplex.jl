# export JuMP_from_SDP_Problem

# see mathprogcomplex : a78593860b662fc1c09dbd3910656b84425729f7

function JuMP_from_SDP_Problem(sdp_pb::MPC.SDP_Problem, mysolver)
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

    JuMP.@constraintref ctrs[1:length(ctr_names)]
    obj = 0
    bodys = Dict()

    # Building constraints and objective bodys
    for ((ctr_name, blockname, var1, var2), f_αβ) in sdp_pb.matrices
        !haskey(bodys, ctr_name) && (bodys[ctr_name] = 0)

        block = sdp_pb.name_to_sdpblock[blockname]
        n = length(block.var_to_id)
        i = block.var_to_id[var1]
        j = block.var_to_id[var2]

        Ai = sparse([i], [j], [1], n, n)
        # if i == j
            # Ai = sparse([i], [j], [1], n, n)
        # else
        #     Ai = sparse([i; j], [j; i], [1; 1], n, n)
        # end

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

    for (i, ctr_name) in enumerate(ctr_names)
        ctrs[i] = JuMP.@constraint(m, bodys[ctr_name] == 0)
    end

    return m

end
