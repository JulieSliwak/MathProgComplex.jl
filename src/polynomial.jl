export Polynomial, print


"""
    Polynomial(poly::SortedDict{Exponent, Number}, degree::Degree)

Define a mathematical polynomial, that is a linear combinaison of product of
`Variable` and conjugated `Variable`.

### Attributes
- `poly` : a dictionary associating an `Exponent` to a a number.
- `degree`: a `Degree` object indicating the global degree of the polynomial :
the maximum of the explicit exponent degrees and of the conjugated exponent
degrees.

### Exemple
```julia
julia > a, b = Variable("a", Complex), Variable("b", Complex)
julia > p = a*b*conj(b) + (2+3im) * b^3
julia > p.poly = SortedDict(a*b*conj(b)=>1, b^3=>2+3im)
julia > p.degree = (3,1)
```
"""
struct Polynomial <: AbstractPolynomial
    poly::SortedDict{Exponent, Number}
    degree::Degree

    function Polynomial(poly::SortedDict{Exponent, Number})
        degexpl, degconj = 0, 0
        for (expo, λ) in poly
            if λ!=0
                degexpl = max(degexpl, expo.degree.explvar)
                degconj = max(degconj, expo.degree.conjvar)
            else
                delete!(poly, expo)
            end
        end
        return new(poly, Degree(degexpl, degconj))
    end
end


#############################
## Constructors
#############################
Polynomial() = Polynomial(SortedDict{Exponent, Number}())


#############################
## Print
#############################
function Base.print(io::IO, P::Polynomial)
  poly = P.poly
  i = length(poly)
  un = Exponent()
  for (expo, λ) in poly

    if expo == un
      print(io, "$λ")
    else
      if λ != 1
        print(io, "(", λ, ")*")
      end
      print(io, expo)
    end
    if i > 1
      print(io, " + ")
    end
    i -= 1
  end
end
