export set_relaxation, get_defaultparams

"""
    relax_ctx = set_relaxation(pb::Problem; ismultiordered=false, issparse=false, symmetries=Set(), hierarchykind=:Complex, renamevars=false, di=Dict{String, Int}(), d=-1)

    Build a `relax_ctx` object containing relaxation choices and problem features : order by constraint, relaxation order by constraint...
"""
function set_relaxation(pb::Problem; ismultiordered::Bool=false,
                                     issparse::Bool=false,
                                     symmetries::Array{DataType, 1}=DataType[],
                                     hierarchykind::Symbol=:Complex,
                                     renamevars::Bool=false,
                                     di::Dict{String, Int}=Dict{String, Int}(),
                                     d::Int=-1,
                                     params=Dict())

    # Check that all variables have a type fitting the hierarchy kind
    for (varname, vartype) in pb.variables
        (vartype<:Int) && error("set_relaxation() : variable $varname,$vartype is integer, unfit for SOS relaxation.\nConsider relaxing it and adding a complementarity constraint.")
        (hierarchykind==:Complex) && !(vartype<:Complex) && error("set_relaxation() : variable $varname,$vartype should be complex for complex hierarchy.")
        (hierarchykind==:Real) && !(vartype<:Real) && error("set_relaxation() : variable $varname,$vartype should be real for real hierarchy.")
    end

    relax_ctx = RelaxationContext()

    relaxparams = relax_ctx.relaxparams

    # Compute each constraint degree
    relctx_setki!(relax_ctx, pb)

    # Store each SDP multiplier type
    relctx_setSDPmulttypes!(relax_ctx, pb, hierarchykind)

    # Relaxation order management
    relctx_setdi!(relax_ctx, pb, di, d)

    rel_ctx_setsymetries!(relax_ctx, pb, symmetries)


    log_POPcharact!(relax_ctx, pb)

    relaxparams[:opt_issparse] = issparse
    relaxparams[:opt_hierarchykind] = hierarchykind
    relaxparams[:opt_multiordered] = ismultiordered
    relaxparams[:opt_globalorder] = relax_ctx.di[get_momentcstrname()]
    (PhaseInvariance in symmetries) && (relaxparams[:opt_sym_phaseinv] = true)

    for (param, val) in params
        if haskey(relax_ctx.relaxparams, param)
            relax_ctx.relaxparams[param] = val
        else
            warn("set_relaxation(): Unhandled parameter $param")
        end
    end

    init_output(relax_ctx)

    print_build_relctx(relax_ctx, pb)
    return relax_ctx
end


function relctx_setki!(relax_ctx, pb::Problem)
    ki = relax_ctx.ki
    ki[get_momentcstrname()] = 0
    for (cstrname, cstr) in pb.constraints
        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)
            ki[cstrname_lo] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
            ki[cstrname_up] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
        else
            ki[get_cstrname(cstrname, cstrtype)] = max(cstr.p.degree.explvar, cstr.p.degree.conjvar)
        end
    end
    return
end

function relctx_setSDPmulttypes!(relax_ctx, pb::Problem, hierarchykind)
    ctrtypes = relax_ctx.cstrtypes
    for (cstrname, cstr) in pb.constraints
        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)
            ctrtypes[cstrname_lo] = (hierarchykind==:Complex ? :SDPC : :SDP)
            ctrtypes[cstrname_up] = (hierarchykind==:Complex ? :SDPC : :SDP)
        elseif cstrtype == :eq
            ctrtypes[get_cstrname(cstrname, cstrtype)] = (hierarchykind==:Complex ? :SymC : :Sym)
        else
            ctrtypes[get_cstrname(cstrname, cstrtype)] = (hierarchykind==:Complex ? :SDPC : :SDP)
        end
    end
    ctrtypes[get_momentcstrname()] = (hierarchykind==:Complex ? :SDPC : :SDP)
    return
end

function relctx_setdi!(relax_ctx, pb::Problem, di, d)
    di_relax = relax_ctx.di
    !((di == Dict{String, Int}()) && (d==-1)) || error("RelaxationContext(): Either di or d should be provided as input.")

    ## Checking order by constraint
    for (cstrname, cstr) in pb.constraints
        cur_order = haskey(di, cstrname) ? di[cstrname] : d

        # Check provided di is suitable wrt constraint degree, add
        cstrtype = get_cstrtype(cstr)
        if cstrtype == :ineqdouble
            cstrname_lo, cstrname_up = get_cstrname(cstrname, cstrtype)
            cur_ki = relax_ctx.ki[cstrname_lo]
            (0 ≤ cur_order-ceil(cur_ki/2)) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value ceil($cur_ki/2).")
            # (cur_ki <= cur_order) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value $cur_ki, hierarchy may be multiordered.")
            di_relax[cstrname_lo] = cur_order
            di_relax[cstrname_up] = cur_order
            # di_relax[cstrname_lo] = max(cur_order, ceil(cur_ki/2))
            # di_relax[cstrname_up] = max(cur_order, ceil(cur_ki/2))
        else # :eq, :ineqlo, :ineqhi
            cur_ki = relax_ctx.ki[get_cstrname(cstrname, cstrtype)]
            (0 ≤ cur_order-ceil(cur_ki/2)) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value ceil($cur_ki/2).")
            # (cur_ki <= cur_order) || warn("RelaxationContext(): Provided order ($cur_order) is lower than constraint $cstrname order ($cur_ki). \nUsing value $cur_ki, hierarchy may be multiordered.")
            di_relax[get_cstrname(cstrname, cstrtype)] = cur_order
            # di_relax[get_cstrname(cstrname, cstrtype)] = max(cur_order, ceil(cur_ki/2))
        end
    end

    ## Setting order for moment contraint(s)
    if haskey(di, get_momentcstrname())
        di_relax[get_momentcstrname()] = di[get_momentcstrname()]
    elseif d!=-1
        di_relax[get_momentcstrname()] = d
    else
        di_relax[get_momentcstrname()] = maximum(values(di_relax))
    end

    # Checking order of objective
    obj_degree = max(pb.objective.degree.explvar, pb.objective.degree.conjvar)
    # if obj_degree > di_relax[get_momentcstrname()]
    #     warn("RelaxationContext(): Moment matrix order $(di_relax[get_momentcstrname()]) is lower than objective degree ($obj_degree). \nUsing value $obj_degree, hierarchy may be multiordered.")
    #     di_relax[get_momentcstrname()] = ceil(obj_degree/2)
    # end

    return
end



function rel_ctx_setsymetries!(relax_ctx, pb, symmetries)
    pbsymmetries = relax_ctx.symmetries
    isa(symmetries, Array) || error("set_relaxation(): symmetries should be an Array of types.")
    for symtype in symmetries
        if has_symmetry(relax_ctx, pb, symtype)
            push!(pbsymmetries, symtype)
        end
    end
    relax_ctx.symmetries = pbsymmetries

    if has_symmetry(relax_ctx, pb, PhaseInvariance)
        relax_ctx.relaxparams[:pb_isphaseinv] = true
    end
    return
end



function log_POPcharact!(relax_ctx::RelaxationContext, pb::Problem)
    params = relax_ctx.relaxparams

    params[:pb_nvar] = length(pb.variables)
    for (varname, vartype) in pb.variables
        if vartype <: Complex
            params[:pb_nvar_cplx] += 1
        elseif vartype <: Real
            params[:pb_nvar_real] += 1
        elseif vartype <: Bool
            params[:pb_nvar_bin] += 1
        end
    end

    params[:pb_nctr] = length(pb.constraints)
    params[:pb_nctr_eq] = -1
    params[:pb_nctr_ineqsimple] = -1
    params[:pb_nctr_ineqdouble] = -1

    params[:pb_maxpolydeg] = maximum(values(relax_ctx.ki))
end

function get_defaultparams()
    defparams = OrderedDict{Symbol, Any}(:pb_name=>"",          # POP parameters
                        :pb_nvar=>0,
                        :pb_nvar_cplx=>0,
                        :pb_nvar_real=>0,
                        :pb_nvar_bin=>0,
                        :pb_nctr=>0,
                        :pb_nctr_eq=>-1,
                        :pb_nctr_ineqsimple=>-1,
                        :pb_nctr_ineqdouble=>-1,
                        :pb_maxpolydeg=>-1,
                        :pb_isphaseinv=>false,
                        :slv_mmtrel_t=>-1.0,                    # Relaxation time and memory consumption
                        :slv_mmtrel_bytes=>-1,
                        :slv_sosrel_t=>-1.0,
                        :slv_sosrel_bytes=>-1,
                        :slv_mskstruct_t=>-1.0,
                        :slv_mskstruct_bytes=>-1,
                        :slv_prosta=>"",                        # SDP solution values
                        :slv_solsta=>"",
                        :slv_primobj=>-1.0,
                        :slv_dualobj=>-1.0,
                        :slv_primfeas=>-1.0,
                        :slv_dualfeas=>-1.0,
                        :slv_solvetime=>-1.0,
                        :opt_hierarchykind=>:Undef,
                        :opt_issparse=>false,                   # Relaxation parameters
                        :opt_multiordered=>false,
                        :opt_globalorder=>-1,
                        :opt_sym_phaseinv=>false,
                        :opt_nb_cliques=>-1,
                        :opt_msk_maxtime=>-1,                   # Default -1 is no time limit
                        :opt_outmode=>0,                        # 0: screen, 1: file, 2: both
                        :opt_outlev=>1,                         # 0: none, 1:summary at moment relaxation, sos relaxation, 2: detailled information, 3: full problems
                        :opt_outname=>"momentsos.log",
                        :opt_outcsv=>0,                         # 0: no csv is written, 1: csv is written
                        :opt_outcsvname=>"momentsos_solve.csv")

    return defparams
end
