
"""
set_constraints!(sdp_pb::SDP_Problem)

Set attributes `obj_keys`, `name_to_ctr`, `id_to_ctr` for a  `sdp_pb`
`SDP_Problem` already set with matrix linear and const contributions.
"""
function set_constraints!(sdp_pb::SDP_Problem)
    # Collecting constraint names
    ctr_names = SortedSet{SDP_CtrObjName}(collect(keys(sdp_pb.cst_ctr)))
    union!(ctr_names, [k[1] for k in keys(sdp_pb.matrices)])
    union!(ctr_names, [k[1] for k in keys(sdp_pb.linear)])

    # Extracting objective keys
    obj_keys = SortedSet{SDP_CtrObjName}()
    for ctr_name in ctr_names
        if (ctr_name[1:2] == ("1", "1"))
            push!(obj_keys, ctr_name)
            delete!(ctr_names, ctr_name)
        end
    end

    length(obj_keys) == 0 && error("No ctrkey matching objective key (\"1\", \"1\", .)")
    sdp_pb.obj_keys = obj_keys

    # Filling constraints attribute
    ctr_id = 1
    for ctr_name in ctr_names
        sdp_pb.name_to_ctr[ctr_name] = (ctr_id, "EQ", 0, 0)
        sdp_pb.id_to_ctr[ctr_id] = ctr_name
        ctr_id += 1
    end
end


"""
set_blocks!(sdp_pb::SDP_Problem)

Set string to id map for each sdp block given a `SDP_Problem` already set with
SDP matrix contributions.
"""
function set_blocks!(sdp_pb::SDP_Problem)
    for (ctr_name, blockname, var1, var2) in keys(sdp_pb.matrices)
        # sanity check
        @assert haskey(sdp_pb.name_to_sdpblock, blockname)

        cur_block = sdp_pb.name_to_sdpblock[blockname]

        # Adding vars and ids to SDP block
        if !haskey(cur_block.var_to_id, var1)
            cur_block.var_to_id[var1] = length(cur_block.var_to_id) + 1
        end
        if !haskey(cur_block.var_to_id, var2)
            cur_block.var_to_id[var2] = length(cur_block.var_to_id) + 1
        end
    end
end

"""
set_scalvars!(sdp_pb::SDP_Problem)

Set string to id map for each scalar variable given a `SDP_Problem` already set
with scalar linear contributions.
"""
function set_scalvars!(sdp_pb::SDP_Problem)
    for (ctr_name, var) in keys(sdp_pb.linear)

        if !haskey(sdp_pb.scalvar_to_id, var)
            sdp_pb.scalvar_to_id[var] = length(sdp_pb.scalvar_to_id) + 1
        end
    end
end