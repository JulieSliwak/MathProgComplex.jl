export build_SDP_Instance_from_SDPPrimal

function build_SDP_Instance_from_SDPPrimal(sdpprimal::SDPPrimal)
    sdp_pb = SDP_Problem()

    # Build numerical descirption of problem, indexed on strings
    set_matrices!(sdp_pb, sdpprimal)
    set_linear!(sdp_pb, sdpprimal)
    set_const!(sdp_pb, sdpprimal)

    # Build structural information and string to id maps
    set_blockvartypes!(sdp_pb, sdpprimal)

    set_constraints!(sdp_pb)
    set_blocks!(sdp_pb)
    set_scalvars!(sdp_pb)

    return sdp_pb
end



function set_matrices!(sdp_pb::SDP_Problem, sdpprimal::SDPPrimal)
    for ((moment, block, γ, δ), coeff) in sdpprimal.blocks
        γ_str, δ_str = format_string(γ), format_string(δ)

        ctr_name = (format_string(moment.conj_part), format_string(moment.expl_part), moment.clique)
        block_name = block
        row_name = max(γ_str, δ_str)
        col_name = min(γ_str, δ_str)

        # sanity check
        @assert !haskey(sdp_pb.matrices, (ctr_name, block_name, row_name, col_name))
        @assert imag(coeff) == 0

        sdp_pb.matrices[(ctr_name, block_name, row_name, col_name)] = coeff
    end
end


function set_linear!(sdp_pb::SDP_Problem, sdpprimal::SDPPrimal)
    # Building symmetric terms into linear terms
    for ((moment, block, var), coeff) in sdpprimal.linsym
        ctr_name = (format_string(moment.conj_part), format_string(moment.expl_part), moment.clique)
        var_name = format_string(var, block)

        # sanity check
        @assert !haskey(sdp_pb.linear, (ctr_name, var_name))
        @assert imag(coeff) == 0

        sdp_pb.linear[(ctr_name, var_name)] = coeff
    end


    # Building scalar linear terms
    for ((moment, var), coeff) in sdpprimal.lin
        ctr_name = (format_string(moment.conj_part), format_string(moment.expl_part), moment.clique)
        var_name = format_string(var)

        # sanity check
        @assert !haskey(sdp_pb.linear, (ctr_name, var_name))
        @assert imag(coeff) == 0

        sdp_pb.linear[(ctr_name, var_name)] = coeff
    end
end

function set_const!(sdp_pb::SDP_Problem, sdpprimal::SDPPrimal)
    for (moment, coeff) in sdpprimal.cst
        ctr_name = (format_string(moment.conj_part), format_string(moment.expl_part), moment.clique)

        # sanity check
        @assert !haskey(sdp_pb.cst_ctr, ctr_name)
        @assert imag(coeff) == 0

        sdp_pb.cst_ctr[ctr_name] = coeff
    end
end


"""
set_blockvartypes!(sdp_pb::SDP_Problem, sdpprimal::SDPPrimal)

Set attributes `name_to_sdpblock` and `id_to_sdpblock` given a primal sdp.
"""
function set_blockvartypes!(sdp_pb::SDP_Problem, sdpprimal::SDPPrimal)
    n_sdp = 0

    for (blockname, blocktype) in sdpprimal.block_to_vartype
        if blocktype == :SDP
            n_sdp += 1

            block = SDP_Block(n_sdp, blockname)
            sdp_pb.name_to_sdpblock[blockname] = block
            sdp_pb.id_to_sdpblock[n_sdp] = block
        elseif blocktype ∉ Set([:Sym])
            error("set_blockvartypes!(): Unknown block variable type '$blocktype' for variable $blockname")
        end
    end
end