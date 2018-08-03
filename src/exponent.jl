export Exponent, show

"""
    Exponent(expo::SortedDict{Variable, Degree})

Define a mathematical exponent, that is a product of `Variable` and conjugated
`Variable`.

### Attributes
- `expo` : a dictionary associating `Variable` objects to a `Degree` objects,
their exponent and the exponent of their conjugates.
- `degree`: a `Degree` object indicating the global degree of the exponent :
the sum of the explicit variables degrees and of the conjugated variables
degrees.

### Exemple
```julia
julia > a, b = Variable("a", Complex), Variable("b", Complex)
julia > Exponent(SortedDict(a=>Degree(1,0), b=>Degree(1,1))) == a*b*conj(b)
```
"""
struct Exponent <: AbstractPolynomial
    expo::SortedDict{Variable, Degree}
    degree::Degree

    function Exponent(expo::SortedDict{Variable, Degree})
        degexpl, degconj = 0,0
        for (var, degree) in expo
            ((degree.explvar < 0) || (degree.conjvar < 0)) && error("Exponent(): Expected non negative exponent for variable $var (got $degree)")
            (isreal(var) && degree.conjvar != 0) && error("Exponent(): Expected nul conj exponent for real variable $var (got $degree)")
            (isbool(var) && degree.explvar ∉ SortedSet([0,1])) && error("Exponent(): Expected boolean exponent for bool $var (got $degree)")
            if degree != Degree(0,0)
                degexpl = max(degexpl, degree.explvar)
                degconj = max(degconj, degree.conjvar)
            else
                delete!(expo, var)
            end
        end
        return new(expo, Degree(degexpl, degconj))
    end
end


#############################
## Constructors
#############################
Exponent() = Exponent(SortedDict{Variable, Degree}())
Exponent(x::Variable) = Exponent(SortedDict{Variable, Degree}(x=>Degree(1,0)))



#############################
## Equality, hash
#############################
function ==(exp1::Exponent, exp2::Exponent)
  (length(exp1.expo) == length(exp2.expo)) || return false
  return isequal(exp1.expo, exp2.expo)
  # for (varname, deg) in exp2.expo
  #   if !haskey(exp1.expo, varname) || exp1.expo[varname] != deg
  #     return false
  #   end
  # end
  # return true
end
!=(exp1::Exponent, exp2::Exponent) = !(exp1 == exp2)
hash(expo::Exponent, h::UInt) = hash(expo.degree, hash(expo.expo, h))



#############################
## Show
#############################
function Base.show(io::IO, exp::Exponent)
  if exp == Exponent()
    print(io, "1")
  else
    expo = exp.expo
    i=length(expo)
    for (var, deg) in expo
      if var.kind <: Complex && deg.conjvar>0
        print(io, "conj(", var, ")")
        if deg.conjvar > 1
          print(io, "^", deg.conjvar)
        end
        if deg.explvar > 0
          print(io, " * ")
        end
      end
      if deg.explvar == 1
        print(io, var)
      elseif deg.explvar > 1
        print(io, var, "^", deg.explvar)
      end
      if i > 1
        print(io, " * ")
      end
      i -= 1
    end
  end
end
