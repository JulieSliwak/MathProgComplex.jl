export build_SDP_Instance_from_SDPDual

function build_SDP_Instance_from_SDPDual(sdpdual::SDPDual)
    sdp_pb = SDP_Problem()

    ## First. Deal fully with moment matrices
    blocktovars = SortedDict{String, SortedSet}()
    n_sdp = 0
    for ((ctrname, cliquename), matrix) in sdpdual.constraints
        if ctrname == get_momentcstrname()
            n_sdp += 1
            blockname = cliquename
            block = SDP_Block(n_sdp, blockname)

            sdp_pb.name_to_sdpblock[blockname] = block
            sdp_pb.id_to_sdpblock[n_sdp] = block

            @assert matrix.matrixkind == :SDP
            blocktovars[blockname] = SortedSet{String}()
            for (α, β) in keys(matrix)
                var1 = format_string(α)
                var2 = format_string(β)

                push!(blocktovars[blockname], var1)
                push!(blocktovars[blockname], var2)
            end
        end
    end
    for (blockname, vars) in blocktovars
        block = sdp_pb.name_to_sdpblock[blockname]
        n=0
        for var in vars
            n += 1
            block.var_to_id[var] = n
        end
    end

    # warn("Done with moment constraints.\n")

    ## Second. Filling obj and ctrs coefficients
    ## Dealing with objective
    ctr_name = ("1", "1", "objective")
    for (moment, fαβ) in sdpdual.objective
        block_name, α, β = split_moment(moment)

        α_str, β_str = format_string(α), format_string(β)
        var1, var2 = max(α_str, β_str), min(α_str, β_str)

        @assert haskey(sdp_pb.name_to_sdpblock, block_name)
        @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var1)
        @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var2)

        if product(α, β) == Exponent()
            !haskey(sdp_pb.cst_ctr, ctr_name) && (sdp_pb.cst_ctr[ctr_name] = 0)

            sdp_pb.cst_ctr[ctr_name] += fαβ
        else
            sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = fαβ * (var1!=var2 ? 0.5 : 1)
        end
    end

    ## Dealing with constraints
    for ((ctrname, cliquename), matrix) in sdpdual.constraints

        if ctrname != get_momentcstrname()

            # Name of auxiliary SDP constrained matrix
            S_name = get_auxSDPmatrix_name(ctrname, cliquename)

            # for all matrix coefficients
            # Deal with S[γ, δ]
            for ((γ, δ), momentpoly) in matrix
                ctr_name = get_dual_ctr_name(ctrname, cliquename, γ, δ)
                γ_str = format_string(γ)
                δ_str = format_string(δ)

                for (moment, fαβ) in momentpoly
                    block_name, α, β = split_moment(moment)

                    α_str, β_str = format_string(α), format_string(β)
                    var1, var2 = max(α_str, β_str), min(α_str, β_str)

                    @assert haskey(sdp_pb.name_to_sdpblock, block_name)
                    @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var1)
                    @assert haskey(sdp_pb.name_to_sdpblock[block_name].var_to_id, var2)

                    # @show (ctr_name, block_name, var1, var2), fαβ
                    if product(α, β) == Exponent()
                        !haskey(sdp_pb.cst_ctr, ctr_name) && (sdp_pb.cst_ctr[ctr_name] = 0)

                        sdp_pb.cst_ctr[ctr_name] += fαβ
                    else
                        sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = fαβ * (var1!=var2 ? 0.5 : 1)
                    end
                end

                # If the matrix is SDP, add link to auxiliary SDP var
                if matrix.matrixkind == :SDP
                    block_name = S_name
                    var1 = max(γ_str, δ_str)
                    var2 = min(γ_str, δ_str)

                    sdp_pb.matrices[(ctr_name, block_name, var1, var2)] = -1 * (var1!=var2 ? 0.5 : 1)
                else
                    # Nothing to do, constraint will be scalar, ==0 by default.
                end
            end
        end
    end

    ## Third. Add auxiliary SDP matrices standing for SDP localizing constraints
    # println("--------------------------------------------------------")
    momentvars = SortedSet(collect(keys(sdp_pb.name_to_sdpblock)))

    for block_name in momentvars
        ctr_name = (block_name*"_1ctr", "1", "1")
        sdp_pb.matrices[(ctr_name, block_name, "1", "1")] = 1
        sdp_pb.cst_ctr[ctr_name] = -1
    end

    # @show momentvars
    # println()
    empty!(blocktovars)

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

            !haskey(blocktovars, blockname) && (blocktovars[blockname] = SortedSet{String}())
            push!(blocktovars[blockname], var1)
            push!(blocktovars[blockname], var2)
        end
    end
    for (blockname, vars) in blocktovars
        block = sdp_pb.name_to_sdpblock[blockname]
        n=0
        for var in vars
            n += 1
            block.var_to_id[var] = n
        end
    end

    ## Deal with equality constaints
    if length(sdpdual.moments_overlap) != 0
        for (expo, cliques) in sdpdual.moments_overlap
            if expo != Exponent()
                @assert length(cliques) > 1

                clique_ref = first(cliques)
                moment_ref = Moment(expo, clique_ref)

                block_name, α, β = split_moment(moment_ref)
                α_str, β_str = format_string(α), format_string(β)
                var1, var2 = max(α_str, β_str), min(α_str, β_str)

                for clique in setdiff(cliques, Set([clique_ref]))
                    ctr_name = get_coupling_ctr_name(clique_ref, clique, var1, var2)

                    sdp_pb.matrices[(ctr_name, clique_ref, var1, var2)] = 1
                    sdp_pb.matrices[(ctr_name, clique, var1, var2)] = -1
                end
            end
        end
    end

    # Build structural information and string to id maps
    # set_blockvartypes!(sdp_pb)

    # set_blocks!(sdp_pb, sdpdual)

    # MathProgComplex.set_blocks!(sdp_pb)

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

function get_coupling_ctr_name(clique_ref, clique, var1, var2)
    return var1, var2, clique_ref*"_"*clique
end


"""
    clique, α, β = split_expo(expo::Exponent)

Split the exponent into two exponents of conjugated and explicit variables in the complex case.
Real case is not supported yet.
"""
function split_moment(moment::Moment)
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

            if α_deg < β_deg
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


function set_blocks!(sdp_pb::SDP_Problem, sdpdual::SDPDual)
    # Collecting moments from each moment constraint
    for ((ctrname, cliquename), matrix) in sdpdual.constraints

        # moment constraint
        if ctrname == get_momentcstrname()
            blockname = cliquename
            block = sdp_pb.name_to_sdpblock[blockname]

            for ((α, β), momentpoly) in matrix
                @assert length(momentpoly) == 1

                moment = first(momentpoly)[1]  #y_αβ
                # @show moment
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
