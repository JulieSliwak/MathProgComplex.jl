# export build_SDP_Instance_from_SDPDual

function build_SDP_Instance_from_SDPDual(sdpdual::MPC.SDPDual)
    sdp_pb = SDP_Problem()

    ## First. Deal fully with moment matrices
    n_sdp = 0
    for ((ctrname, cliquename), matrix) in sdpdual.constraints
        if ctrname == get_momentcstrname()
            n_sdp += 1
            blockname = cliquename
            block = SDP_Block(n_sdp, blockname)

            sdp_pb.name_to_sdpblock[blockname] = block
            sdp_pb.id_to_sdpblock[n_sdp] = block

            @assert matrix.matrixkind == :SDP
            for (α, β) in keys(matrix)
                var1 = format_string(α)
                if !haskey(block.var_to_id, var1)
                    block.var_to_id[var1] = length(block.var_to_id)+1
                end

                var2 = format_string(β)
                if !haskey(block.var_to_id, var2)
                    block.var_to_id[var2] = length(block.var_to_id)+1
                end
            end

            println(" * $blockname:\n$block\n")
        end
    end

    warn("Done with moment constraints.\n")

    ## Second. Filling obj and ctrs coefficients
    ## Dealing with objective
    ctr_name = ("1", "1", "objective")
    for (moment, fαβ) in sdpdual.objective
        block_name, α, β = split_moment(moment)

        α_str, β_str = format_string(α), format_string(β)
        var1, var2 = min(α_str, β_str), max(α_str, β_str)

        @assert haskey(sdp_pb.name_to_sdpblock, block_name)
        @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var1)
        @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var2)
        @show (ctr_name, block_name, var1, var2), fαβ
        sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = -fαβ
    end

    warn("Done with moment objective.\n")

    ## Dealing with constraints
    for ((ctrname, cliquename), matrix) in sdpdual.constraints
        if ctrname != get_momentcstrname()

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
                    block_name, α, β = split_moment(moment)

                    α_str, β_str = format_string(α), format_string(β)
                    var1, var2 = min(α_str, β_str), max(α_str, β_str)

                    @assert haskey(sdp_pb.name_to_sdpblock, block_name)
                    @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var1)
                    @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var2)

                    @show (ctr_name, block_name, var1, var2), fαβ
                    sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = fαβ
                end

                # If the matrix is SDP, add link to auxiliary SDP var
                if matrix.matrixkind == :SDP
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

    ## Third. Add auxiliary SDP matrices standing for SDP localizing constraints
    println("--------------------------------------------------------")
    momentvars = SortedSet(collect(keys(sdp_pb.name_to_sdpblock)))

    @show momentvars
    println()

    n_sdp = length(sdp_pb.name_to_sdpblock)
    for ((ctr_name, blockname, var1, var2), f_αβ) in sdp_pb.matrices
        if blockname ∉ momentvars
            if !haskey(sdp_pb.name_to_sdpblock, blockname)
                n_sdp += 1
                block = SDP_Block(n_sdp, blockname)

                sdp_pb.name_to_sdpblock[blockname] = block
                sdp_pb.id_to_sdpblock[n_sdp] = block
            else
                block = sdp_pb.name_to_sdpblock[blockname]
            end

            @assert !haskey(block.var_to_id, var1)
            @assert !haskey(block.var_to_id, var2)
            if !haskey(block.var_to_id, var1)
                block.var_to_id[var1] = length(block.var_to_id)+1
            end

            if !haskey(block.var_to_id, var2)
                block.var_to_id[var2] = length(block.var_to_id)+1
            end

            println(" * $blockname:\n$block\n")
        end
    end

    ## Deal with equality constaints
    if length(sdpdual.moments_overlap) != 0
        error("Coupling constraints not implemented yet.")
    end

    # Build structural information and string to id maps
    # set_blockvartypes!(sdp_pb)

    # set_blocks!(sdp_pb, sdpdual)

    # MathProgComplex.set_blocks!(sdp_pb)


    # warn("222222222222222222222222222")
    # println(sdp_pb)
    # warn("222222222222222222222222222")

    MathProgComplex.set_constraints!(sdp_pb)
    MathProgComplex.set_scalvars!(sdp_pb)

    return sdp_pb
end


function get_dual_ctr_name(ctrname, cliquename, γ, δ)
    return replace(ctrname, " ", "")*"_"*replace(cliquename, " ", ""), format_string(γ), format_string(δ)
end

function get_auxSDPmatrix_name(ctrname, cliquename)
    return "S_"*replace(ctrname, " ", "")*"_"*replace(cliquename, " ", "")
end

"""
    α, β = split_expo(expo::Exponent)

    Split the exponent into two exponents of conjugated and explicit variables in the complex case.
    Real case is not supported yet.
"""
function split_moment(moment::MPC.Moment)
    α, β = Exponent(), Exponent()

    if moment.conj_part != Exponent()
        # Complex moment, separation is trivial
        α = moment.conj_part
        β = moment.expl_part

    else
        # Real moment, more complicated...
        explsum, conjsum = MathProgComplex.get_sumdegs(moment.expl_part)
        @assert conjsum == 0

        α_deg = β_deg = 0
        α_curadd = β_curadd = 0
        for (var, degree) in moment.expl_part
            @assert degree.conjvar == 0

            if α_curadd < β_curadd
                α_curadd = ceil(degree.explvar / 2)
                β_curadd = floor(degree.explvar / 2)
            else
                α_curadd = floor(degree.explvar / 2)
                β_curadd = ceil(degree.explvar / 2)
            end

            product!(α, Exponent(SortedDict(var=>Degree(α_curadd, 0))))
            α_deg += α_curadd
            product!(β, Exponent(SortedDict(var=>Degree(β_curadd, 0))))
            β_deg += β_curadd
        end

        @assert product(α, β) == product(moment.expl_part, moment.conj_part)
        αexplsum, αconjsum = get_sumdegs(α)
        βexplsum, βconjsum = get_sumdegs(β)

        @assert min(αexplsum, βexplsum) == floor(explsum/2)
        @assert max(αexplsum, βexplsum) == ceil(explsum/2)
    end


    return moment.clique, α, β
end

"""
set_blockvartypes!(sdp_pb::SDP_Problem)

# Set attributes `name_to_sdpblock` and `id_to_sdpblock` given a primal sdp.
"""
function set_blockvartypes!(sdp_pb::MPC.SDP_Problem)
    n_sdp = 0

    blocknames = SortedSet([k[2] for k in keys(sdp_pb.matrices)])

    for blockname in blocknames
        n_sdp += 1

        block = SDP_Block(n_sdp, blockname)
        sdp_pb.name_to_sdpblock[blockname] = block
        sdp_pb.id_to_sdpblock[n_sdp] = block
    end
end


function set_blocks!(sdp_pb::MPC.SDP_Problem, sdpdual::MPC.SDPDual)
    # Collecting moments from each moment constraint
    for ((ctrname, cliquename), matrix) in sdpdual.constraints

        # moment constraint
        if ctrname == get_momentcstrname()
            blockname = cliquename
            block = sdp_pb.name_to_sdpblock[blockname]

            for ((α, β), momentpoly) in matrix
                @assert length(momentpoly) == 1

                moment = first(momentpoly)[1]  #y_αβ
                @show moment
                α_str = format_string(moment.conj_part)
                β_str = format_string(moment.expl_part)

                if !haskey(block.var_to_id, α_str)
                    block.var_to_id[α_str] = length(block.var_to_id) + 1
                end
                if !haskey(block.var_to_id, β_str)
                    block.var_to_id[β_str] = length(block.var_to_id) + 1
                end
            end

        end
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