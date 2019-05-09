export ==, !=, hash, isless

#############################
## Degree
#############################
function isless(deg1::Degree, deg2::Degree)
  return (deg1.explvar<deg2.explvar) || (deg1.explvar==deg2.explvar && deg1.conjvar<deg2.conjvar)
end


#############################
## Variable
#############################
function isless(a::Variable, b::Variable)
  return isless(a.name, b.name)
end


#############################
## Exponent
#############################
function isless(exp1::Exponent, exp2::Exponent)
    ## Handle one exponent, empty structures are not handled by following general sort alg.
    if length(exp1) == 0        # if exp1 == 1
        return length(exp2) > 0
    end

    if length(exp2) == 0        # exp1!=1, if exp2 == 1
        return false
    end

    iter_result1 = iterate(exp1)
    iter_result2 = iterate(exp2)
    while iter_result1 !== nothing  && iter_result2 !== nothing
          (i1, state1) = iter_result1
          (i2, state2) = iter_result2
        if isequal(i1, i2)
            # continue
        else
            return isless(i1,i2)
        end
        iter_result1 = iterate(exp1, state1)
        iter_result2 = iterate(exp2, state2)

        if iter_result1 === nothing && iter_result2 === nothing
            return false
        elseif iter_result1 === nothing && iter_result2 !== nothing
            return true
        elseif iter_result1 !== nothing && iter_result2 === nothing
            return false
        end
    end
    return false
end
"""
    isless_degree(exp1::Exponent, exp2::Exponent)

    Order sorting elements on their sum of degrees at a first level.
    Test show performance is comparable to previous sorting function.
"""
function isless_degree(exp1::Exponent, exp2::Exponent)
    # First order level: sum of degrees
    exp1_deg = 0
    iter_result1 = iterate(exp1)
    while iter_result1 !== nothing
        (i1, state1) = iter_result1
        exp1_deg += i1[2].explvar + i1[2].conjvar
        iter_result1 = iterate(exp1, state1)
    end

    exp2_deg = 0
    iter_result2 = iterate(exp2)
    while iter_result2 !== nothing
        (i2, state2) = iter_result2
        exp2_deg += i2[2].explvar + i2[2].conjvar
        iter_result2 = iterate(exp2, state2)
    end

    if exp1_deg < exp2_deg
        return true
    elseif exp1_deg > exp2_deg
        return false
    end

    # Second order level: sort with variables and degrees
    iter_result1 = iterate(exp1)
    iter_result2 = iterate(exp2)
    while iter_result1 !== nothing && iter_result2 !== nothing
        (i1, state1) = iter_result1
        (i2, state2) = iter_result2

        if isequal(i1, i2)
            #continue
        else
            return isless(i1,i2)
        end
        iter_result1 = iterate(exp1, state1)
        iter_result2 = iterate(exp2, state2)

        if iter_result1 === nothing && iter_result2 === nothing
            return false
        elseif iter_result1 === nothing && iter_result2 !== nothing
            return true
        elseif iter_result1 !== nothing && iter_result2 === nothing
            return false
        end
    end
    return false
end

#############################
## Polynomial
#############################
# NOTE: order on polynomials ?
