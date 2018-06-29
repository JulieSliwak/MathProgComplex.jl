"""
    momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem)

    Build the sparsitty pattern and variables decomposition for laying out the moment or SOS hierarchy
"""
function build_sparsity(relax_ctx::RelaxationContext, problem::Problem, max_cliques::Dict{String, Set{Variable}})

    ((relax_ctx.issparse == false) && (length(max_cliques) > 1)) && error("build_sparsity(): Relaxation is not sparse, one clique is expected (not $(length(max_cliques)))")

    # Build localizing constraints order and variable set.
    localizingmat_param = Dict{String, Tuple{Set{String}, Int}}()
    for (ctrname, ctr) in problem.constraints
        ctrtype = get_cstrtype(ctr)
        ctrcliques = get_locctrcliques(ctr.p, max_cliques)

        if ctrtype == :ineqdouble
            ctrname_lo, ctrname_up = get_cstrname(ctrname, ctrtype)
            di_lo, ki_lo = relax_ctx.di[ctrname_lo], relax_ctx.ki[ctrname_lo]
            di_up, ki_up = relax_ctx.di[ctrname_up], relax_ctx.ki[ctrname_up]

            if relax_ctx.hierarchykind==:Complex
                localizingmat_param[ctrname_lo] = (ctrcliques, di_lo-ki_lo)
                localizingmat_param[ctrname_up] = (ctrcliques, di_up-ki_up)
            elseif relax_ctx.hierarchykind==:Real
                localizingmat_param[ctrname_lo] = (ctrcliques, di_lo-ceil(ki_lo/2))
                localizingmat_param[ctrname_up] = (ctrcliques, di_up-ceil(ki_up/2))
            else
                error("build_sparsity(): Unknown relaxation kind $(relax_ctx.hierarchykind). Should be `:Real` or `:Complex`")
            end
        else # :ineqlo, :ineqhi, :eq
            di, ki = relax_ctx.di[get_cstrname(ctrname, ctrtype)], relax_ctx.ki[get_cstrname(ctrname, ctrtype)]
            if relax_ctx.hierarchykind == :Complex
                localizingmat_param[get_cstrname(ctrname, ctrtype)] = (ctrcliques, di-ki)
            elseif relax_ctx.hierarchykind == :Real
                localizingmat_param[get_cstrname(ctrname, ctrtype)] = (ctrcliques, di-ceil(ki/2))
            else
                error("build_sparsity(): Unknown relaxation kind $(relax_ctx.hierarchykind). Should be `:Real` or `:Complex`")
            end
        end
    end

    # Build moment constraints order and variable set.
    momentmat_param = Dict{String, Int}()
    for (cliquename, cliquevars) in max_cliques
        cur_d::Int = -1
        for (ctrname, (ctrcliques, _)) in localizingmat_param
            if length(ctrcliques) == 1 && cliquename == first(ctrcliques)
                cur_d = max(cur_d, relax_ctx.di[ctrname])
            end
        end
        if cur_d == -1 # Clique does not match a full constraint
            cur_d = minimum(values(relax_ctx.di))
        end
        momentmat_param[cliquename] = cur_d
    end

    return momentmat_param, localizingmat_param
end


"""
pvars = get_variables(p)

Collect all variables appearing in `p`.
"""
function get_variables(p::Polynomial)
    pvars = Set{Variable}()
    for (expo, coeff) in p
        for (var, deg) in expo
            push!(pvars, var)
        end
    end
    return pvars
end

"""
    locctrcliques = get_locctrcliques(p, max_cliques)

    Find a minimal set of cliques gathering all variables from polynomial `p`.
"""
function get_locctrcliques(p::Polynomial, max_cliques::Dict{String, Set{Variable}})
    ctrvars = get_variables(p)

    # Build constraint variables to cliques dict
    var_to_cliques = Dict{Variable, Set{String}}()
    for (clique, cliquevars) in max_cliques
        for var in intersect(cliquevars, ctrvars)
            if !haskey(var_to_cliques, var)
                var_to_cliques[var] = Set{String}()
            end
            push!(var_to_cliques[var], clique)
        end
    end

    def_cliques = Set{String}()
    unaffected_vars = Set{Variable}(keys(var_to_cliques))

    keepon = true
    i = 0
    while keepon
        i += 1
        # Find which variables appear in one clique only, remove them from unaffected_vars
        for var in unaffected_vars
            cliques = var_to_cliques[var]
            if length(cliques) == 1
                push!(def_cliques, first(cliques))
                delete!(unaffected_vars, var)
            end
        end

        # Remove variables involved in at least one clique previously selected
        for var in unaffected_vars
            inter = intersect(var_to_cliques[var], def_cliques)
            if !isempty(inter)
                var_to_cliques[var] = Set{String}([first(inter)]) # NOTE: proper way to choose in inter here ?
                delete!(unaffected_vars, var)
            end
        end

        # Hopefully all variables are treated that way. Else repeat this process by choosing a clique. Again, which one ?
        if length(unaffected_vars) != 0
            # warn("get_locctrcliques(): length(unaffected_vars) = $(length(unaffected_vars))") # TODO: better logging system...
            cliques_from_unaffvar = Dict{String, Int}()
            for var in unaffected_vars
                for clique in var_to_cliques[var]
                    # haskey(cliques_from_unaffvar, clique) || (cliques_from_unaffvar[clique] = 0)
                    # cliques_from_unaffvar[clique] += 1
                    addindex!(cliques_from_unaffvar, 1, clique)
                end
            end
            cur_clique, cur_pop = first(cliques_from_unaffvar)[1], first(cliques_from_unaffvar)[2]
            for (clique, pop) in cliques_from_unaffvar
                if pop > cur_pop
                    cur_clique = clique
                    cur_pop = pop
                end
            end
            push!(def_cliques, cur_clique)
        else
            keepon = false
        end
    end

    return def_cliques
end

function get_maxcliques(relax_ctx, problem)
    vars = Set{Variable}([Variable(name, kind) for (name, kind) in problem.variables])
    return Dict{String, Set{Variable}}("clique1"=>vars)
end

function get_WB5cliques(relax_ctx, problem)
    if !relax_ctx.issparse
        return get_maxcliques(relax_ctx, problem)
    else
        maxcliques = Dict{String, Set{Variable}}()
        maxcliques["clique1"] = Set{Variable}([
            Variable("BaseCase_1_VOLT_Im", Real),
            Variable("BaseCase_1_VOLT_Re", Real),
            Variable("BaseCase_2_VOLT_Im", Real),
            Variable("BaseCase_2_VOLT_Re", Real),
            Variable("BaseCase_3_VOLT_Im", Real),
            Variable("BaseCase_3_VOLT_Re", Real)])
        maxcliques["clique2"] = Set{Variable}([
            Variable("BaseCase_2_VOLT_Im", Real),
            Variable("BaseCase_2_VOLT_Re", Real),
            Variable("BaseCase_3_VOLT_Im", Real),
            Variable("BaseCase_3_VOLT_Re", Real),
            Variable("BaseCase_4_VOLT_Im", Real),
            Variable("BaseCase_4_VOLT_Re", Real),
            Variable("BaseCase_5_VOLT_Im", Real),
            Variable("BaseCase_5_VOLT_Re", Real)])
        return maxcliques
    end
end

function get_case9cliques(relax_ctx, problem)
    if !relax_ctx.issparse
        return get_maxcliques(relax_ctx, problem)
    else
        maxcliques = Dict{String, Set{Variable}}()
        maxcliques["clique1"] = Set{Variable}([
            Variable("BaseCase_1_VOLT_Im", Real),
            Variable("BaseCase_1_VOLT_Re", Real),
            Variable("BaseCase_5_VOLT_Im", Real),
            Variable("BaseCase_5_VOLT_Re", Real),
            Variable("BaseCase_4_VOLT_Im", Real),
            Variable("BaseCase_4_VOLT_Re", Real),
            Variable("BaseCase_9_VOLT_Im", Real),
            Variable("BaseCase_9_VOLT_Re", Real),
            Variable("BaseCase_8_VOLT_Im", Real),
            Variable("BaseCase_8_VOLT_Re", Real)])
        maxcliques["clique2"] = Set{Variable}([
            Variable("BaseCase_2_VOLT_Im", Real),
            Variable("BaseCase_2_VOLT_Re", Real),
            Variable("BaseCase_9_VOLT_Im", Real),
            Variable("BaseCase_9_VOLT_Re", Real),
            Variable("BaseCase_8_VOLT_Im", Real),
            Variable("BaseCase_8_VOLT_Re", Real),
            Variable("BaseCase_7_VOLT_Im", Real),
            Variable("BaseCase_7_VOLT_Re", Real),
            Variable("BaseCase_6_VOLT_Im", Real),
            Variable("BaseCase_6_VOLT_Re", Real)])
        maxcliques["clique3"] = Set{Variable}([
            Variable("BaseCase_3_VOLT_Im", Real),
            Variable("BaseCase_3_VOLT_Re", Real),
            Variable("BaseCase_7_VOLT_Im", Real),
            Variable("BaseCase_7_VOLT_Re", Real),
            Variable("BaseCase_6_VOLT_Im", Real),
            Variable("BaseCase_6_VOLT_Re", Real),
            Variable("BaseCase_5_VOLT_Im", Real),
            Variable("BaseCase_5_VOLT_Re", Real),
            Variable("BaseCase_4_VOLT_Im", Real),
            Variable("BaseCase_4_VOLT_Re", Real)])
        return maxcliques
    end
end



"""
    vars, blocname = collect_cliquesvars(clique_keys, max_cliques)

    Collect variables of `cliques_keys` cliques, described in `max_cliques`
"""
function collect_cliquesvars(clique_keys::Set{String}, max_cliques::Dict{String, Set{Variable}})
    # Collect variables involved in constraint
    vars = Set{Variable}()
    blocname = ""
    for clique_key in clique_keys
        union!(vars, max_cliques[clique_key])
        blocname = blocname*clique_key*"_"
    end
    return vars, blocname[1:end-1]
end

function print(io::IO, max_cliques::Dict{String, Set{Variable}})
    for cliquename in sort(collect(keys(max_cliques)))
        vars = max_cliques[cliquename]

        print(io, "$cliquename = ")
        for var in vars print(io, "$var, ") end
        @printf(io, "\b\b \n")
    end
end


#################################################################################
## Old stuff
"""
    SparsityPattern

    Type for storing and working on sparsitty patterns.
"""
type SparsityPattern end

"""
    sparsity_pattern = compute_sparsitypattern(problem, di, ki)
    Compute the sparsity_pattern corresponding to the given partial orders.
"""
function compute_sparsitypattern(problem::Problem, relax_ctx)
    println("\n=== compute_sparsitypattern(problem::Problem, relax_ctx)")
    println("Compute the sparsity_pattern corresponding to the given partial orders.")
    println("-> Nb vertex / edges:      xx / xx")
    println("-> Order by vertex:        xx / xx (mean/std)")
    sparsity_pattern = SparsityPattern()
    return sparsity_pattern::SparsityPattern
end

"""
    compute_chordalextension!(sparsity_pattern)
    Compute the chordal extension on the provided sparsity_pattern.
"""
function compute_chordalextension!(sparsity_pattern::SparsityPattern)
    println("\n=== compute_chordalextension!(sparsity_pattern::SparsityPattern)")
    println("Compute the chordal extension on the provided sparsity_pattern.")
    println("-> Nb added edges:                 xx")
    println("-> Nb vertex / edges:              xx / xx")
    println("-> Order by vertex:                xx / xx (mean/std)")
end

"""
    maxcliques = compute_maxcliques(sparsity_pattern)
    Compute a `Array{Set{Variable}}` describing the maximum cliques on the provided sparsity_pattern.
"""
function compute_maxcliques(sparsity_pattern::SparsityPattern)
    println("\n=== compute_maxcliques(sparsity_pattern::SparsityPattern)")
    println("Compute a `Array{Set{Variable}}` describing the maximum cliques on the provided sparsity_pattern.")
    println("-> Nb cliques / nb vertex:         xx / xx")
    println("-> Nb vertex by clique:            xx / xx (mean/std)")
    maxcliques = Array{Set{Variable}}()
    return maxcliques::Array{Set{Variable}}
end


"""
    varsbycstr = compute_varsbycstr(problem)
    Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in each constraint.
"""
function compute_varsbycstr(problem::Problem)
    println("\n=== compute_varsbycstr(problem::Problem)")
    println("Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in each constraint.")
    println("-> Nb poly constraints:                xx")
    println("-> Nb variables by poly constraint:    xx / xx (mean/std)")
    varsbycstr = Dict{String, Set{Variable}}()
    return varsbycstr::Dict{String, Set{Variable}}
end

"""
    cliquevarsbycstr = compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.
"""
function compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)
    println("\n=== compute_varsbycstr(sparsity_pattern, max_cliques, varsbycstr)")
    println("Compute a `Dict{String, Set{Variable}}` providing the set of variables involved in the SDP localizing matrix corresponding to each constraint.")
    println("-> Nb SDP constraints:                 xx")
    println("-> Nb variables by SDP constraint:     xx / xx (mean/std)")
    cliquevarsbycstr = Dict{String, Set{Variable}}()
    return cliquevarsbycstr::Dict{String, Set{Variable}}
end

"""
    orderbyclique = compute_cliqueorders(sparsity_pattern, di, varsbycstr, max_cliques)
    Compute a `Array{Int}` providing the relaxation order corresponding to each clique.
"""
function compute_cliqueorders(sparsity_pattern, varsbycstr, max_cliques, relax_ctx)
    println("\n=== compute_cliqueorders(sparsity_pattern, varsbycstr, max_cliques, relax_ctx)")
    println("Compute a `Array{Int}` providing the relaxation order corresponding to each clique.")
    println("-> Nb cliques:                         xx")
    println("-> Relaxation order by clique:         xx / xx (mean/std)")
    orderbyclique = Dict{Int, Int}()
    return orderbyclique::Dict{Int, Int}
end
