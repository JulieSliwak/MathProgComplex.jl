# """
#     mm = MomentMatrix(vars::SortedSet{Variable}, d, symmetries)

#     Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
#     Only monomials featuring all `symmetries` appear in the moment matrix.
# """
function MomentMatrix{T}(relax_ctx::RelaxationContext, vars::Set{Variable},
                                                       d::Int,
                                                       symmetries::Set{DataType},
                                                       matrixkind::Symbol;
                                                       default_clique::String="",
                                                       var_to_cliques::Dict{Variable, Set{String}}=Dict{Variable, Set{String}}()) where T<:Number

    mm = Dict{Tuple{Exponent, Exponent}, Dict{Moment, T}}()

    @assert default_clique!="" || var_to_cliques!=Dict{Variable, Set{String}}()

    ## Computing exponents for available variables
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    for cexp in conjexpos
        for rexp in realexpos
            expo = product(cexp, rexp)

            ## Checking if current exponent has required symmetries
            issym = true
            for sym in symmetries
                issym = issym && has_symmetry(relax_ctx, expo, sym)
            end

            ## Storing only lower triangular matrix
            if issym && cexp ≥ rexp

                # Get exponent clique
                expo_clique = default_clique
                if var_to_cliques!=Dict{Variable, Set{String}}()
                    expo_clique = get_exponentclique(expo, var_to_cliques)
                end

                mm[(cexp, rexp)] = Dict{Moment, T}(Moment(expo, expo_clique)=>convert(T, 1))
            end
        end
    end
    return MomentMatrix{T}(mm, Set(vars), d, matrixkind)
end

# function copy(mm::MomentMatrix)
#     return MomentMatrix(copy(mm.mm), mm.vars, mm.order, mm.matrixkind)
# end

function print(io::IO, mm::MomentMatrix{T}) where T<:Number
    keylen = maximum(x->length("($(x[1]), $(x[2])) → "), keys(mm.mm))

    maxorder = 0

    for key in sort(collect(keys(mm.mm)))
        momentpoly = mm.mm[key]

        print_string(io, "($(key[1]), $(key[2])) → ", keylen)

        momentpoly_sortedkeys = sort(collect(keys(momentpoly)))
        moment_first = pop!(momentpoly_sortedkeys)
        val_first = momentpoly[moment_first]
        maxorder = max(maxorder, moment_first.expl_part.degree.explvar)

        println(io, "$(moment_first.clique) -- $(moment_first.conj_part) × $(moment_first.expl_part) × $val_first")

        for moment in momentpoly_sortedkeys
            val = momentpoly[moment]
            maxorder = max(maxorder, moment.expl_part.degree.explvar)

            (moment_first, val_first) == (moment, val) && continue
            println(io, " "^(keylen+1), "$(moment.clique) -- $(moment.conj_part) × $(moment.expl_part) × $val")
        end
    end
    info(io, "maxorder is $maxorder")
    print(io, " $(mm.matrixkind)")
end



"""
    cliquename = get_exponentclique(expo, var_to_cliques)

    Determine which clique expo fits in, that is which cliques contain all variables of expo.
    Error if no such clique are found.
"""
function get_exponentclique(expo::Exponent, var_to_cliques::Dict{Variable, Set{String}})
    cliques = Set{String}()

    ## If expo is one, return default clique
    expo == Exponent() && return "clique_un"

    union!(cliques, var_to_cliques[first(expo)[1]])
    for (var, deg) in expo
        cliques = intersect(cliques, var_to_cliques[var])
    end

    length(cliques) == 0 && error("get_exponentclique(): $expo is split amongst several cliques.\nMaximal cliques provided are not suitable for this relaxation.")

    clique = first(cliques)
    # length(cliques) > 1 && warn("get_exponentclique(): $expo appears in $(length(cliques)) cliques : $cliques.\nChoosing first one $clique") ## TODO: better logging system
    return clique
end

# ##########################
# ## Moment matrix algebra
# ##########################

## AbstractPolynomial types
function product!(mm::MomentMatrix{M}, p::T, var_to_cliques::Dict{Variable, Set{String}}) where T<:Union{AbstractPolynomial, Number} where M<:Number
    for (key, momentpoly) in mm.mm
        mm.mm[key] = product(momentpoly, p, var_to_cliques)
    end
    return nothing
end

function product(momentpoly::Dict{Moment, M}, p::T, var_to_cliques::Dict{Variable, Set{String}}) where T<:Union{AbstractPolynomial, Number} where M<:Number
    return product(momentpoly, convert(Polynomial, p), var_to_cliques)
end

function product(momentpoly::Dict{Moment, M}, p::Polynomial, var_to_cliques::Dict{Variable, Set{String}}) where M<:Number
    resmpoly = Dict{Moment, M}()

    for (expo, val1) in p
        for (moment, val2) in momentpoly
            resmoment = product(moment, expo, var_to_cliques)

            # haskey(resmpoly, resmoment) || (resmpoly[resmoment] = convert(M, 0.0))
            # resmpoly[resmoment] += val1*val2
            addindex!(resmpoly, val1*val2, resmoment)

            isnull(resmpoly[resmoment]) && delete!(resmpoly, resmoment)
        end
    end

    return resmpoly
end

function product(moment::Moment, expo::Exponent, var_to_cliques::Dict{Variable, Set{String}})
    resexpo = product(moment.conj_part, moment.expl_part)
    product!(resexpo, expo)

    clique = get_exponentclique(resexpo, var_to_cliques)
    return Moment(resexpo, clique)
end

# # function evaluate(mm::MomentMatrix, pt::Point)
# #     mm_eval = SortedDict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
# #     for (key, p) in mm.mm
# #         res = evaluate(p, pt)
# #         if res == Polynomial()
# #             delete!(mm_eval, key)
# #         else
# #             mm_eval[key] = res
# #         end
# #     end
# #     return MomentMatrix(mm_eval, setdiff(mm.vars, SortedSet(keys(pt))), mm.order, mm.matrixkind)
# # end
