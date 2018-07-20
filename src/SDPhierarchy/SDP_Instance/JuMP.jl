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

        # cstr1 = JuMP.@constraint(m, [i=1:n,j=1:n], var[i, j] == var[i+n, j+n])
        # # for i in 1:n, j in 1:n
        # #     println("$i, $j  --> $(cstr1[i, j])")
        # # end
        # cstr2 = JuMP.@constraint(m, [i=1:n,j=1:i], var[n+i, j] == - var[n+j, i])
        # # for i in 1:n, j in 1:i
        # #     println("$i, $j --> $(cstr2[i, j])")
        # # end
    end

    varlin = Dict()
    for varname in keys(sdp_pb.scalvar_to_id)
        varlin[varname] = JuMP.@variable(m, basename=varname)
    end

    ## Constraints
    # Constraint storage
    ctr_names = SortedSet([k[1] for k in keys(sdp_pb.matrices)])
    union!(ctr_names, [k[1] for k in keys(sdp_pb.linear)])
    union!(ctr_names, [k for k in keys(sdp_pb.cst_ctr)])
    ctr_names = collect(ctr_names)

    warn("$(length(ctr_names)) constraints found.")
    JuMP.@constraintref ctrs[1:length(ctr_names)]

    bodys = Dict()

    for ((ctr_name, blockname, var1, var2), f_αβ) in sdp_pb.matrices
        !haskey(bodys, ctr_name) && (bodys[ctr_name] = 0)

        block = sdp_pb.name_to_sdpblock[blockname]
        n = length(block.var_to_id)
        i = block.var_to_id[var1]
        j = block.var_to_id[var2]

        Ai = sparse([i], [j], 1, n, n)
        bodys[ctr_name] += vecdot( Zi[blockname], Ai) * f_αβ
    end

    for ((ctr_name, varname), f_αβ) in sdp_pb.linear
        !haskey(bodys, ctr_name) && (bodys[ctr_name] = 0)

        bodys[ctr_name] += varlin[varname] * f_αβ
    end

    for (ctr_name, f_αβ) in sdp_pb.cst_ctr
        bodys[ctr_name] += f_αβ
    end

    # Laying out constraints
    for (i, ctr_name) in enumerate(ctr_names)
        ctrs[i] = JuMP.@constraint(m, bodys[ctr_name] == 0)
    end

    return m

end
