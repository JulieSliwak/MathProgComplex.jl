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
    # exp1_deg = sum(get_sumdegs(exp1))
    # exp2_deg = sum(get_sumdegs(exp2))

    exp1_deg = 0
    state1 = start(exp1)
    while !done(exp1, state1)
        (i1, state1) = next(exp1, state1)
        exp1_deg += i1[2].explvar + i1[2].conjvar
    end

    exp2_deg = 0
    state2 = start(exp2)
    while !done(exp2, state2)
        (i2, state2) = next(exp2, state2)
        exp2_deg += i2[2].explvar + i2[2].conjvar
    end

    # First order level: sum of degrees
    if exp1_deg < exp2_deg
        return true
    elseif exp1_deg > exp2_deg
        return false
    end

    state1 = start(exp1)
    state2 = start(exp2)

    while !done(exp1, state1) && !done(exp2, state2)
        (i1, state1) = next(exp1, state1)
        (i2, state2) = next(exp2, state2)

        if isequal(i1, i2)
            continue
        else
            return isless(i1,i2)
        end
        if done(exp1, state1) && done(exp2, state2)
            return false
        elseif done(exp1, state1) && !done(exp2, state2)
            return true
        elseif !done(exp1, state1) && done(exp2, state2)
            return false
        end
    end

    return false
end

#############################
## Polynomial
#############################
# NOTE: order on polynomials ?
