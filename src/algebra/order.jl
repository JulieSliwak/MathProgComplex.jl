export ==, !=, hash, isless

#############################
## Degree
#############################
==(d1::Degree, d2::Degree) = (d1.explvar==d2.explvar) && (d1.conjvar==d2.conjvar)
!=(d1::Degree, d2::Degree) = !(d1 == d2)

hash(deg::Degree, h::UInt) = hash(deg.explvar, hash(deg.conjvar, h))

function isless(deg1::Degree, deg2::Degree)
  return (deg1.explvar<deg2.explvar) || (deg1.explvar==deg2.explvar && deg1.conjvar<deg2.conjvar)
end


#############################
## Variable
#############################
==(x::Variable, y::Variable) = (x.name==y.name) && (x.kind==y.kind)
!=(x::Variable, y::Variable) = !(x == y)

hash(var::Variable, h::UInt) = hash(var.name, hash(var.kind, h))

function isless(a::Variable, b::Variable)
  return isless(a.name, b.name)
end



#############################
## Exponent
#############################
function ==(exp1::Exponent, exp2::Exponent)
  (length(exp1.expo) == length(exp2.expo)) || return false

  for (varname, deg) in exp2.expo
    if !haskey(exp1.expo, varname) || exp1.expo[varname] != deg
      return false
    end
  end
  return true
end
!=(exp1::Exponent, exp2::Exponent) = !(exp1 == exp2)

hash(expo::Exponent, h::UInt) = hash(expo.degree, hash(expo.expo, h))

"""
  isless(exp1, exp2)

  BEWARE: Order is valid for fully real or conjugated exponents. Mixed case is not treated.
"""
function isless(exp1::Exponent, exp2::Exponent)
  exp1_explsum, exp1_conjsum = get_sumdegs(exp1)
  exp2_explsum, exp2_conjsum = get_sumdegs(exp2)
  # exp1_explsum==0 || exp1_conjsum==0 || warn("isless(::Exponent, Exponent): exp1 has expl and conj vars, order may be ill defined...") # TODO
  # exp2_explsum==0 || exp2_conjsum==0 || warn("isless(::Exponent, Exponent): exp2 has expl and conj vars, order may be ill defined...")
  exp1_deg = exp1_explsum + exp1_conjsum
  exp2_deg = exp2_explsum + exp2_conjsum
  if exp1_deg < exp2_deg
    return true
  elseif exp1_deg == exp2_deg
    vars = SortedSet(keys(exp1))
    union!(vars, keys(exp2))
    for var in vars
      if !haskey(exp2, var) # hence haskey(exp1, var)
        return true
      elseif !haskey(exp1, var) # hence haskey(exp2, var)
        return false
      elseif exp1[var] > exp2[var]
        return true
      elseif exp1[var] == exp2[var]
        continue
      else
        return false
      end
    end
  else
    return false
  end
  return false
end



#############################
## Polynomial
#############################
function ==(p1::Polynomial, p2::Polynomial)
  if (length(p1) != length(p2)) || (p1.degree != p2.degree)
    return false
  end
  for (exponent, coeff) in p1
    if !haskey(p2, exponent) || (coeff != p2.poly[exponent])
      return false
    end
  end
  true
end
!=(pol1::Polynomial, pol2::Polynomial) = !(pol1 == pol2)

# NOTE: order on polynomials ?



#############################
## Point
#############################
function ==(pt1::Point, pt2::Point)
  if length(pt1) != length(pt2)
    return false
  end
  for (var1, val1) in pt1
    if !haskey(pt2, var1) || pt2[var1] != val1
      return false
    end
  end
  true
end
!=(pt1::Point, pt2::Point) = !(pt1 == pt2)
