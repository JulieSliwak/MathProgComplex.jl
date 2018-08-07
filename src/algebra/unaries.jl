export convert, conj, real, imag, abs2, norm

## convert functions
function convert(::Type{Polynomial}, expo::Exponent)
    return Polynomial(SortedDict{Exponent, Number}(expo=>1.0))
end

function convert(::Type{Polynomial}, x::Variable)
    return Polynomial(SortedDict{Exponent, Number}(Exponent(SortedDict(x=>Degree(1,0)))=>1.0))
end

function convert(::Type{Polynomial}, λ::Number)
    return Polynomial(SortedDict{Exponent, Number}(Exponent()=>λ))
end

function convert(::Type{Point}, pt::SortedDict{Variable, T}) where T
    return Point(convert(SortedDict{Variable, Number}, pt))
end

function convert(::Type{AbstractPolynomial}, λ::Number)
    return Polynomial(SortedDict{Exponent, Number}(Exponent()=>λ))
end


## Conjugate
function conj(d::Degree)
    return Degree(d.conjvar, d.explvar)
end

function conj(x::Variable)
    if x.kind<:Complex
        return Exponent(SortedDict(x=>Degree(0,1)))
    else
        return Exponent(SortedDict(x=>Degree(1,0)))
    end
end

function conj(expo::Exponent)
    expodict = SortedDict{Variable, Degree}()
    for (var, deg) in expo.expo
        if iscomplex(var)
            expodict[var] = conj(deg)
        else
            expodict[var] = deg
        end
    end
    return Exponent(expodict)
end

function conj(p::Polynomial)
    pdict = SortedDict{Exponent, Number}()
    for (expo, λ) in p
        pdict[conj(expo)] = conj(λ)
    end
    return Polynomial(pdict)
end


## Real, imag
real(p::Polynomial) = (p + conj(p))/2
imag(p::Polynomial) = (p - conj(p))/(2im)

function real(p::T) where T<:AbstractPolynomial
    return (p + conj(p))/2
end

function imag(p::T) where T<:AbstractPolynomial
    return (p + -conj(p))/(2im)
end


## Norm, abs2
function abs2(p::T) where T<:AbstractPolynomial
    return p*conj(p)
end

function norm(pt::Point, p::Real = 2)
    p > 0 || error(LOGGER, "norm(): p should be positive ($p ≤ 0).")
    if p == Inf
        maxi = -Inf
        for (var, val) in pt
            # iscomplex(var) && warn("norm(): var $(var) is not real, splitting real and imag parts.")
            if !iscomplex(var) && (val > maxi)
                maxi = val
            elseif iscomplex(var) && max(real(val), imag(val)) > maxi
                maxi = max(real(val), imag(val))
            end
        end
        return maxi
    elseif p == 1
      s = 0
      for (var, val) in pt
          if !iscomplex(var)
            s += abs(val)
          else
            # warn("norm(): var $(var) is not real, splitting real and imag parts.")
            s += abs(real(val)) + abs(imag(val))
          end
      end
      return s
    else
      s = 0
      for (var, val) in pt
          if !iscomplex(var)
            s += val^p
          else
            # warn("norm(): var $(var) is not real, splitting real and imag parts.")
            s += real(val)^p + imag(val)^p
          end
      end
      return s^(1/p)
  end
end