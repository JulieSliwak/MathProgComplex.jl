function build_SDP_Instance_from_SDPDual(sdpdual::SDPDual)
    sdp_pb = SDP_Problem()

    ## Dealing with objective
    ctr_name = ("1", "1", "objective")
    for (moment, fαβ) in sdpdual.objective
        block_name = moment.clique
        α = format_string(moment.conj_part)
        β = format_string(moment.expl_part)

        var1, var2 = min(α, β), max(α, β)
        # @show (ctr_name, block_name, var1, var2), fαβ
        sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = -fαβ
    end

    ## Dealing with constraints
    for ((ctrname, cliquename), matrix) in sdpdual.constraints
        if ctrname == get_momentcstrname()
            # Moment constraint, do nothing !
        else
            # Localizing constraint

            # Name of auxiliary SDP constrained matrix
            S_name = get_auxSDPmatrix_name(ctrname, cliquename)
            info("Sname: ", S_name)

            # for all matrix coefficients
            # Deal with S[γ, δ]
            for ((γ, δ), momentpoly) in matrix
                ctr_name = get_dual_ctr_name(ctrname, cliquename, γ, δ)
                γ_str = format_string(γ)
                δ_str = format_string(δ)

                for (moment, fαβ) in momentpoly
                    block_name = moment.clique
                    α = format_string(moment.conj_part)
                    β = format_string(moment.expl_part)

                    var1, var2 = min(α, β), max(α, β)
                    # @show (ctr_name, block_name, var1, var2), fαβ
                    sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = fαβ
                end

                # If the matrix is SDP, add link to auxiliary SDP var
                if matrix.matrixkind == :SDP
                    warn("SDP, ctr")
                    block_name = S_name
                    var1 = min(γ_str, δ_str)
                    var2 = max(γ_str, δ_str)

                    # @show (ctr_name, block_name, var1, var2), -1
                    sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = -1
                else
                    # Nothing to do, constraint will be scalar, ==0 by default.
                end
            end
        end
    end

    ## Deal with equality constaints
    if length(sdpdual.moments_overlap) != 0
        error("Coupling constraints not implemented yet.")
    end

    # Build structural information and string to id maps
    set_blockvartypes!(sdp_pb)

    MathProgComplex.set_constraints!(sdp_pb)
    MathProgComplex.set_blocks!(sdp_pb)
    MathProgComplex.set_scalvars!(sdp_pb)

    return sdp_pb
end


function get_dual_ctr_name(ctrname, cliquename, γ, δ)
    return replace(ctrname, " ", "")*"_"*replace(cliquename, " ", ""), format_string(γ), format_string(δ)
end

function get_auxSDPmatrix_name(ctrname, cliquename)
    return replace(ctrname, " ", "")*"_"*replace(cliquename, " ", "")
end


"""
set_blockvartypes!(sdp_pb::SDP_Problem)

# Set attributes `name_to_sdpblock` and `id_to_sdpblock` given a primal sdp.
"""
function set_blockvartypes!(sdp_pb::SDP_Problem)
    n_sdp = 0

    blocknames = SortedSet([k[2] for k in keys(sdp_pb.matrices)])

    for blockname in blocknames
        n_sdp += 1

        block = SDP_Block(n_sdp, blockname)
        sdp_pb.name_to_sdpblock[blockname] = block
        sdp_pb.id_to_sdpblock[n_sdp] = block
    end
end



# mutable struct MomentMatrix{T}
#     mm::DictType{Tuple{Exponent, Exponent}, DictType{Moment, T}}
#     vars::Set{Variable}
#     order::Int
#     matrixkind::Symbol            # Either :SDP or :Sym
# end

# include(joinpath("base_types", "momentmatrix.jl"))

# """
#     momentrel = SDPDual(obj, cstrs, moment_overlap)

#     Store a Moment Relaxation problem.
# """
# struct SDPDual{T}
#     objective::DictType{Moment, T}                                  # A linear comb. of moments, to be maximized
#     constraints::DictType{Tuple{String, String}, MomentMatrix{T}}   # A set of moment matrices, either SDP or Null. A constraint (`key[1]`) can be split on several cliques (`key[2]`)
#     moments_overlap::DictType{Exponent, Set{String}}                # A set of clique per exponent, describing coupling constraints
# end