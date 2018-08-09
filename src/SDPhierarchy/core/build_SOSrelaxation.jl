export build_SOSrelaxation

function build_SOSrelaxation(relaxctx::RelaxationContext, mmtrelax_pb::SDPDual{T}; debug=false) where T<:Number
    sdpblocks = DictType{Tuple{CtrName, String, Exponent, Exponent}, T}()
    sdplinsym = DictType{Tuple{CtrName, String, Exponent}, T}()
    sdplin = DictType{Tuple{CtrName, Exponent}, T}()
    sdpcst = DictType{CtrName, T}()
    block_to_vartype = DictType{String, Symbol}()

    ## Build blocks dict
    for ((cstrname, cliquename), mmt) in mmtrelax_pb.constraints
        block_name = get_blockname(cstrname, cliquename, mmtrelax_pb)
        if mmt.matrixkind == :Null
            block_to_vartype[block_name] = (relaxctx.hierarchykind == :Real)?:Sym:SymC
        else
            block_to_vartype[block_name] = mmt.matrixkind
        end

        for ((γ, δ), poly) in mmt.mm
            for (moment, λ) in poly
                debug && (expo = product(moment.conj_part, moment.expl_part))
                debug && (@assert expo.degree.explvar == moment.expl_part.degree.explvar)
                debug && (@assert expo.degree.conjvar == moment.conj_part.degree.conjvar)
                # Check the current monomial has correct degree
                if (relaxctx.hierarchykind==:Complex) && ((moment.expl_part.degree.explvar > relaxctx.di[cstrname]) || (moment.conj_part.degree.conjvar > relaxctx.di[cstrname]))
                    warn(LOGGER, "build_SOSrelaxation(): Found exponent pair of degree $(expo.degree) > $(relaxctx.di[cstrname]) for Complex hierarchy.\n($(expo), at $((γ, δ)) of MM matrix)")
                elseif (relaxctx.hierarchykind==:Real) && ((moment.expl_part.degree.explvar > 2*relaxctx.di[cstrname]) || (moment.conj_part.degree.conjvar != 0))
                    warn(LOGGER, "build_SOSrelaxation(): Found exponent pair of degree $(expo.degree) > 2*$(relaxctx.di[cstrname]) for Real hierarchy.\n($(expo), at $((γ, δ)) of MM matrix)")
                end
                !isnan(λ) || warn(LOGGER, "build_SOSrelaxation(): isNaN ! constraint $cstrname - clique $blocname - mm entry $((γ, δ)) - moment $(moment)")

                # Add the current coeff to the SDP problem
                # Constraints are fα - ∑ Bi.Zi = 0
                if mmt.matrixkind == :SDP || mmt.matrixkind == :SDPC
                    key = (moment, block_name, γ, δ)
                    debug && (@assert !haskey(sdpblocks, key))

                    sdpblocks[key] = -λ
                elseif mmt.matrixkind == :Null
                    key = (moment, block_name, product(γ, δ))
                    val = -λ * (γ!=δ ? 2 : 1)

                    haskey(sdplinsym, key) || (sdplinsym[key] = 0)
                    sdplinsym[key] += val
                else
                    error(LOGGER, "build_SOSrelaxation(): Unhandled matrix kind $(mmt.matrixkind) for ($cstrname, $cliquename)")
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
        if !haskey(sdpcst, moment)
            sdpcst[moment] = 0.0
        end

        sdpcst[moment] += fαβ
    end

    sosrelaxation = SDPPrimal{T}(block_to_vartype, sdpblocks, sdplinsym, sdplin, sdpcst)

    print_build_SOSrelax(relaxctx, sosrelaxation)
    return sosrelaxation
end
