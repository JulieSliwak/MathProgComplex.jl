"""
    momentrelaxation = MomentRelaxation(relax_ctx, problem, moment_param::Dict{String, Tuple{Set{String}, Int}}, max_cliques::Dict{String, Set{Variable}})

    Compute the `momentrelaxation` of `problem` corresponding to the clique decomposition `max_cliques` and parameters `moment_param`.
"""
function MomentRelaxation{T}(relax_ctx::RelaxationContext, problem::Problem,
                                                           momentmat_param::Dict{String, Int},
                                                           localizingmat_param::Dict{String, Tuple{Set{String}, Int}},
                                                           max_cliques::Dict{String, Set{Variable}}) where T<:Number

    var_to_cliques = Dict{Variable, Set{String}}()
    for (clique, vars) in max_cliques
        for var in vars
            haskey(var_to_cliques, var) || (var_to_cliques[var] = Set{String}())
            push!(var_to_cliques[var], clique)
        end
    end

    ## Building linear-in-moments objective
    objective = Dict{Moment, T}()
    for (expo, val) in problem.objective
        clique = get_exponentclique(expo, var_to_cliques)
        objective[Moment(expo, clique)] = val
    end


    ## Building linear matrix inequalities
    momentmatrices = Dict{Tuple{String, String}, MomentMatrix{T}}()

    ## Build moment matrix
    for (cliquename, vars) in max_cliques
        dcl = momentmat_param[cliquename]
        momentmatrices[(get_momentcstrname(), cliquename)] = MomentMatrix{T}(relax_ctx, vars, dcl, relax_ctx.symmetries,
                                                                                                   relax_ctx.cstrtypes[get_momentcstrname()],
                                                                                                   default_clique = cliquename)
    end

    ## Build localizing matrices
    for (cstrname, cstr) in problem.constraints

        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname_lower(cstrname), get_cstrname_upper(cstrname)

            # Deal with lower inequality
            clique_keys, order = localizingmat_param[cstrname_lo]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            # length(clique_keys) == 1 || error("MomentRelaxation(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")

            mmt = MomentMatrix{T}(relax_ctx, vars, order, relax_ctx.symmetries,
                                                          relax_ctx.cstrtypes[cstrname_lo],
                                                          var_to_cliques = var_to_cliques)
            # print_with_color(:green, "$cstrname, :Lo\n") ##NOTE: find better logging system.
            product!(mmt, cstr.p - cstr.lb, var_to_cliques)
            momentmatrices[(cstrname_lo, cliquename)] = mmt

            # Deal with upper inequality
            clique_keys, order = localizingmat_param[cstrname_up]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            # length(clique_keys) == 1 || error("MomentRelaxation(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")

            mmt = MomentMatrix{T}(relax_ctx, vars, order, relax_ctx.symmetries,
                                                          relax_ctx.cstrtypes[cstrname_up],
                                                          var_to_cliques = var_to_cliques)
            # print_with_color(:green, "$cstrname, :Up\n") ##NOTE: find better logging system.
            product!(mmt, cstr.ub - cstr.p, var_to_cliques)
            momentmatrices[(cstrname_up, cliquename)] = mmt

            # # Deal with upper inequality, no recomputing of variables or moment matrix if possible
            # clique_keys_up, order_up = localizingmat_param[cstrname_up]
            # if collect(clique_keys) != collect(clique_keys_up)
            #     warn("clique keys different from lower and upper side of double constraint")
            #     length(clique_keys_up) == 1 || error("MomentRelaxation(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")
            #     vars, cliquename = collect_cliquesvars(clique_keys_up, max_cliques)

            #     mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries,
            #                                                   relax_ctx.cstrtypes[cstrname_hi],
            #                                                   var_to_cliques = var_to_cliques)
            # elseif order_up != order
            #     warn("order different from lower and upper side of double constraint")
            #     mmt = MomentMatrix(relax_ctx, vars, order_up, relax_ctx.symmetries,
            #                                                   relax_ctx.cstrtypes[cstrname_hi],
            #                                                   var_to_cliques = var_to_cliques)
            # end


        else
            # either cstrtype == :ineqlo, :ineqhi, :eq
            clique_keys, order = localizingmat_param[get_cstrname(cstrname, cstrtype)]
            vars, cliquename = collect_cliquesvars(clique_keys, max_cliques)
            # length(clique_keys) == 1 || error("MomentRelaxation(): constraint $cstrname spans several cliques ($clique_keys).\nNot supported yet.")

            mmt = MomentMatrix{T}(relax_ctx, vars, order, relax_ctx.symmetries,
                                                          relax_ctx.cstrtypes[get_cstrname(cstrname, cstrtype)],
                                                          var_to_cliques = var_to_cliques)


            # print_with_color(:green, "$cstrname, :Up\n") ##NOTE: find better logging system.
            product!(mmt, get_normalizedpoly(cstr, cstrtype), var_to_cliques)

            momentmatrices[(get_cstrname(cstrname, cstrtype), cliquename)] = mmt
        end
    end

    ## Locate clique overlapping moments
    expo_to_cliques = Dict{Exponent, Set{String}}()

    # Collect Exponents per clique (moment matrix)
    for ((ctrobj, clique), mmtmat) in momentmatrices
        if ctrobj == get_momentcstrname()

            for (key, momentpoly) in mmtmat.mm
                for (moment, coeff) in momentpoly
                    expo = product(moment.conj_part, moment.expl_part)
                    haskey(expo_to_cliques, expo) || (expo_to_cliques[expo] = Set{String}())
                    push!(expo_to_cliques[expo], moment.clique)
                end
            end
        end
    end

    nb_expos = length(expo_to_cliques)
    for (expo, cliques) in expo_to_cliques
        length(cliques) > 1 || delete!(expo_to_cliques, expo)
    end

    momentrelaxation = MomentRelaxation{T}(objective, momentmatrices, expo_to_cliques)

    print_build_momentrelax(relax_ctx, momentrelaxation, nb_expos)

    return momentrelaxation
end
