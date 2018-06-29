"""
    hassym = has_symmetry(relax_ctx, pb::Problem, symtype::T) where T

    Check the `pb` problem for general symmetry `sym`, by examining the objective polynomial and then all constraints bodys.
"""
function has_symmetry(relax_ctx, pb::Problem, sym::T) where T
    # Check objective function for symmetry sym
    if !has_symmetry(relax_ctx, pb.objective, sym)
        return false
    end

    # Check constraints for symmetry sym
    for (cstrname, cstr) in pb.constraints
        if !has_symmetry(relax_ctx, cstr.p, sym)
            return false
        end
    end
    return true
end

"""
    hassym = has_symmetry(relax_ctx, p::Polynomial, sym::T) where T

    Check the `p` polynomial for general symmetry `sym`, by checking each exponent for `sym` symmetry.
"""
function has_symmetry(relax_ctx, p::Polynomial, sym::T) where T
    for (expo, λ) in p
        if !has_symmetry(relax_ctx, expo, sym)
            return false
        end
    end
    return true
end



## Phase Invariance symmetry
"""
    has_symmetry(relax_ctx, expo::Exponent, sym::T) where T<:Type{PhaseInvariance}

    Check whether `expo` is phase invariant, i.e. the exponent value is the same when evaluated at a complex point or a phase shifted complex point (*-1 for real variables).
"""
function has_symmetry(relax_ctx, expo::Exponent, sym::T) where T<:Type{PhaseInvariance}
    explsum, conjsum = get_sumdegs(expo)
    if relax_ctx.hierarchykind == :Real
        (conjsum != 0) && error("is_homogeneous(): Exponent $expo has conjugated variables for a real hierarchy.")
        return explsum % 2 == 0
    elseif relax_ctx.hierarchykind == :Complex
        return explsum == conjsum
    else
        error("is_homogeneous(expo, kind): kind should be either :Real or :Complex ($kind here).")
    end
end