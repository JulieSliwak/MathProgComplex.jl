get_cstrname_lower(cstrname::String) = cstrname*"_lo"
get_cstrname_upper(cstrname::String) = cstrname*"_hi"
get_cstrname_eq(cstrname::String) = cstrname*"_eq"

get_momentcstrname() = "moment_cstr"

function get_blockname(cstrname, cliquename, mmtrelax_pb)
    cstrcliques = Set(Iterators.filter(x->x[1] == cstrname, keys(mmtrelax_pb.constraints)))
    if length(cstrcliques) == 1
        return cstrname
    else
        cstrname == get_momentcstrname() || warn("get_blockname(): several cliques found for $cstrname : $cstrcliques")
        return cstrname*"_"*cliquename
    end
end

function get_cstrname(cstrname::String, cstrtype::Symbol)
    if cstrtype == :ineqdouble
        return get_cstrname_lower(cstrname), get_cstrname_upper(cstrname)
    elseif cstrtype == :eq
        return get_cstrname_eq(cstrname)
    elseif cstrtype == :ineqlo
        return cstrname
        # return get_cstrname_lower(cstrname)
    elseif cstrtype == :ineqhi
        return cstrname
        # return get_cstrname_upper(cstrname)
    else
        error("get_cstrname(): Unknown type $cstrtype.")
    end
end

function get_normalizedpoly(cstr::Constraint, cstrtype::Symbol)
    if cstrtype == :eq
        return cstr.p - cstr.lb
    elseif cstrtype == :ineqlo
        return cstr.p - cstr.lb
    elseif cstrtype == :ineqhi
        return cstr.ub - cstr.p
    else
        error("get_normalizedpoly(): Unhandeld type $cstrtype.")
    end
end

function get_pbcstrname(cstrname::String)
    if ismatch(r".+(\_lo|\_hi|\_eq)", cstrname)
        return cstrname[1:end-3]
    else
        return cstrname
    end
end

function get_ccmultvar(relaxctx::RelaxationContext, moment::Exponent, clique1::String, clique2::String)
    if relaxctx.hierarchykind == :Real
        return Variable("lagmult_cc_$(format_string(moment))_$(clique1)_$(clique2)", Real)
    elseif relaxctx.hierarchykind == :Complex
        return Variable("lagmult_cc_$(format_string(moment))_$(clique1)_$(clique2)", Complex)
    else
        error("get_ccmultvar(): Unhandled hierarchy kind $(relaxctx.hierarchykind) ")
    end
end

function format_string(α::Exponent, β::Exponent)
    s = "$α,$β"
    return replace(s, " ", "")
end

function format_string(s1::String, s2::String)
    s = "$s1,$s2"
    return replace(s, " ", "")
end

function format_string(α::Exponent)
    s = "$α"
    return replace(s, " ", "")
end

function format_string(α::Exponent, block::String)
    s = "$(α)_$block"
    return format_string(α)*"_$block"
end

shortname_moment(n::Int) = "mmt_$n"

function change_eq_to_ineq!(problem::Problem)
    for (ctrname, ctr) in problem.constraints
        if get_cstrtype(ctr) == :eq
            rm_constraint!(problem, ctrname)
            add_constraint!(problem, get_cstrname_lower(ctrname), 0 << (ctr.p - ctr.lb))
            add_constraint!(problem, get_cstrname_upper(ctrname), 0 << (ctr.ub - ctr.p))
        end
    end
end


## Misc
function print_cmat(mat::AbstractArray, round = 1e-3)
    for i=1:size(mat, 1)
        for j=1:size(mat, 2)
            re, im = real(mat[i, j]), imag(mat[i, j])
            @printf("% 5.4f", re)
            @printf(" ")
            @printf("%+5.4fim", im)
            @printf("   ")
        end
        @printf("\n")
    end
end

import Base: hashindex, isslotempty, isslotmissing, isslotfilled, _setindex!, rehash!



### Efficient dict manipulation
"""
    addindex!(h::Dict{K,V}, v0, key::K) where V<:Number where K

    Add the value `v0` to `h[key]` if `key` is already a key, else add the pair `key=>v0`.
"""
function addindex!(h::Dict{K,V}, v0, key::K) where V<:Number where K
    v = convert(V, v0)
    index = ht_keyindex2!(h, key)

    if index > 0
        h.age += 1
        @inbounds h.keys[index] = key
        @inbounds h.vals[index] += v
    else
        @inbounds _setindex!(h, v, key, -index)
    end

    return h
end




#####################################################################################
###  Utils functions, duplicate from base/dict.jl ... (TODO: fix this ducplicate)
#####################################################################################

# get the index where a key is stored, or -pos if not present
# and the key would be inserted at pos
# This version is for use by setindex! and get!
function ht_keyindex2!(h::Dict{K,V}, key) where V where K
    maxallowedprobe = 16
    maxprobeshift   = 6

    age0 = h.age
    sz = length(h.keys)
    iter = 0
    maxprobe = h.maxprobe
    index = hashindex(key, sz)
    avail = 0
    keys = h.keys

    @inbounds while true
        if isslotempty(h,index)
            if avail < 0
                return avail
            end
            return -index
        end

        if isslotmissing(h,index)
            if avail == 0
                # found an available slot, but need to keep scanning
                # in case "key" already exists in a later collided slot.
                avail = -index
            end
        elseif key === keys[index] || isequal(key, keys[index])
            return index
        end

        index = (index & (sz-1)) + 1
        iter += 1
        iter > maxprobe && break
    end

    avail < 0 && return avail

    maxallowed = max(maxallowedprobe, sz>>maxprobeshift)
    # Check if key is not present, may need to keep searching to find slot
    @inbounds while iter < maxallowed
        if !isslotfilled(h,index)
            h.maxprobe = iter
            return -index
        end
        index = (index & (sz-1)) + 1
        iter += 1
    end

    rehash!(h, h.count > 64000 ? sz*2 : sz*4)

    return ht_keyindex2!(h, key)
end