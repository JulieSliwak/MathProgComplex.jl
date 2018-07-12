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
"""
  isless(exp1, exp2)

  BEWARE: Order is valid for fully real or conjugated exponents. Mixed case is not treated.
"""
function isless(exp1::Exponent, exp2::Exponent)
    # return isless(exp1.expo, exp2.expo)
  exp1_explsum, exp1_conjsum = get_sumdegs(exp1)
  exp2_explsum, exp2_conjsum = get_sumdegs(exp2)
  # exp1_explsum==0 || exp1_conjsum==0 || warn("isless(::Exponent, Exponent): exp1 has expl and conj vars, order may be ill defined...") # TODO
  # exp2_explsum==0 || exp2_conjsum==0 || warn("isless(::Exponent, Exponent): exp2 has expl and conj vars, order may be ill defined...")
  exp1_deg = exp1_explsum + exp1_conjsum
  exp2_deg = exp2_explsum + exp2_conjsum
  if exp1_deg < exp2_deg
    return true
  elseif exp1_deg == exp2_deg
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
  #   vars = SortedSet(keys(exp1))
  #   union!(vars, keys(exp2))
  #   for var in vars
  #     if !haskey(exp2, var) # hence haskey(exp1, var)
  #       return true
  #     elseif !haskey(exp1, var) # hence haskey(exp2, var)
  #       return false
  #     elseif exp1[var] > exp2[var]
  #       return true
  #     elseif exp1[var] == exp2[var]
  #       continue
  #     else
  #       return false
  #     end
  #   end
  # else
  #   return false
  end
  return false
end


#############################
## Polynomial
#############################
# NOTE: order on polynomials ?
