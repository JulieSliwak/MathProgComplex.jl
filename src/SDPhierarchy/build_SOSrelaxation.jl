function build_SOSrelaxation(relaxctx::RelaxationContext, mmtrelax_pb::MomentRelaxation{T}; debug=false) where T<:Number
    sdpblocks = Dict{Tuple{Moment, String, Exponent, Exponent}, T}()
    sdplinsym = Dict{Tuple{Moment, String, Exponent}, T}()
    sdplin = Dict{Tuple{Moment, Exponent}, T}()
    sdpcst = Dict{Moment, T}()
    block_to_vartype = Dict{String, Symbol}()

    ## Build blocks dict
    for ((cstrname, cliquename), mmt) in mmtrelax_pb.constraints
        block_name = get_blockname(cstrname, cliquename, mmtrelax_pb)
        block_to_vartype[block_name] = mmt.matrixkind

        for ((γ, δ), poly) in mmt.mm
            for (moment, λ) in poly
                debug && (expo = product(moment.conj_part, moment.expl_part))
                debug && (@assert expo.degree.explvar == moment.expl_part.degree.explvar)
                debug && (@assert expo.degree.conjvar == moment.conj_part.degree.conjvar)
                # Check the current monomial has correct degree
                if (relaxctx.hierarchykind==:Complex) && ((moment.expl_part.degree.explvar > relaxctx.di[cstrname]) || (moment.conj_part.degree.conjvar > relaxctx.di[cstrname]))
                    warn("build_SOSrelaxation(): Found exponent pair of degree $(expo.degree) > $(relaxctx.di[cstrname]) for Complex hierarchy.\n($(expo), at $((γ, δ)) of MM matrix)")
                elseif (relaxctx.hierarchykind==:Real) && ((moment.expl_part.degree.explvar > 2*relaxctx.di[cstrname]) || (moment.conj_part.degree.conjvar != 0))
                    warn("build_SOSrelaxation(): Found exponent pair of degree $(expo.degree) > 2*$(relaxctx.di[cstrname]) for Real hierarchy.\n($(expo), at $((γ, δ)) of MM matrix)")
                end
                !isnan(λ) || warn("build_SOSrelaxation(): isNaN ! constraint $cstrname - clique $blocname - mm entry $((γ, δ)) - moment $(moment)")

                # Add the current coeff to the SDP problem
                # Constraints are fα - ∑ Bi.Zi = 0
                if mmt.matrixkind == :SDP || mmt.matrixkind == :SDPC
                    key = (moment, block_name, γ, δ)
                    debug && (@assert !haskey(sdpblocks, key))

                    sdpblocks[key] = -λ
                elseif mmt.matrixkind == :Sym || mmt.matrixkind == :SymC
                    ## TODO: look into ht_keyindex2!, ht_keyindex for avoiding two dict table lookup
                    # Maybe implement this operation (if haskey add, else set) using 'setindex!(h::Dict{K,V}, v0, key::K) where V where K' as inspiration
                    key = (moment, block_name, product(γ, δ))
                    val = -λ * (γ!=δ ? 2 : 1)

                    addindex!(sdplinsym, val, key)

                    # haskey(sdplinsym, key) || (sdplinsym[key] = 0)
                    # sdplinsym[key] += val
                else
                    error("build_SOSrelaxation(): Unhandled matrix kind $(mmt.matrixkind) for ($cstrname, $cliquename)")
                end

            end
        end
    end

    ## Build linear dict
    # Enforce clique coupling constraints on moments
    for (expo, cliques) in mmtrelax_pb.moments_overlap
        expo == Exponent() && continue
        # print_with_color(:light_cyan, "-> $expo - $(collect(cliques))\n")
        @assert length(cliques)>1
        cliqueref = first(cliques)

        refmoment = Moment(expo, cliqueref)
        for clique in setdiff(cliques, [cliqueref])
            curmoment = Moment(expo, clique)
            var = Exponent(get_ccmultvar(relaxctx, expo, cliqueref, clique))

            @assert !haskey(sdplin, (refmoment, var))
            @assert !haskey(sdplin, (curmoment, var))
            sdplin[(refmoment, var)] =  1
            sdplin[(curmoment, var)] = -1
        end
    end

    ## Build constants dict
    for (moment, fαβ) in mmtrelax_pb.objective
        # Determine which moment to affect the current coefficient.

        # Constraints are fα - ∑ Bi.Zi = 0
        addindex!(sdpcst, fαβ, moment)

        # if !haskey(sdpcst, moment)
        #     sdpcst[moment] = 0.0
        # end

        # # Constraints are fα - ∑ Bi.Zi = 0
        # sdpcst[moment] += fαβ
    end

    sosrelaxation = SDPInstance{T}(block_to_vartype, sdpblocks, sdplinsym, sdplin, sdpcst)

    print_build_SOSrelax(relaxctx, sosrelaxation)
    return sosrelaxation
end


# """
#     α, β = split_expo(expo::Exponent)

#     Split the exponent into two exponents of conjugated and explicit variables in the complex case.
#     Real case is not supported yet.
# """
# function split_expo(relaxctx::RelaxationContext, expo::Exponent)
#     α, β = Exponent(), Exponent()

#     for (var, deg) in expo
#         product!(α, Exponent(Dict(var=>Degree(0, deg.conjvar))))
#         product!(β, Exponent(Dict(var=>Degree(deg.explvar, 0))))
#     end

#     if (relaxctx.hierarchykind == :Real) && (α.degree != Degree(0,0))
#         error("split_expo(): Inconsistent degree $α, $β found for $(relaxctx.hierarchykind) hierarchy.")
#     end
#     return α, β
# end
