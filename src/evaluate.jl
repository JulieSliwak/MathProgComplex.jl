export evaluate

# NOTE: Point are sparse vectors, hence non referenced variables are assumed to be null.
function evaluate(p::Polynomial, pt::Point; partial=false)
    res=0
    expldeg = conjdeg = 0
    for (expo, λ) in p
        res += λ*evaluate(expo, pt; partial=partial)
    end
    typeof(res)<:Polynomial && update_degree!(res)

    # Fully evaluated polynomial
    if typeof(res) == Polynomial && res.degree == Degree(0,0)
        return first(res)[2]
    end
    return res
end

function evaluate(expo::Exponent, pt::Point; partial=false)
    res=1
    for (var, deg) in expo.expo
        res *= (evaluate(var, pt; partial=partial)^deg.explvar) * (conj(evaluate(var, pt; partial=partial))^deg.conjvar)
    end
    return res
end

function evaluate(x::Variable, pt::Point; partial=false)
    if !haskey(pt, x)
        return partial ? x : 0.0
    elseif x.kind<:Complex
        return pt[x]
    else
        return real(pt[x])
    end
end
