export Constraint, <<, >>, ==

"""
    Constraint(p::AbstractPolynomial, lb::Number, ub::Number)

Define a mathematical complex constraint, from a polynomial body and two complex
bounds. Equivalent to the two real inequality produced from real and imag part
of `p`, `lb` and `ub`.

### Attributes
- `p::Polynomial` : the polynomial body of the constraint.
- `lb::Number` : complex lower bound.
- `ub::Number` : complex upper bound.
- `precond::Symbol` : a symbol describing the preconditioning method to be
applied for this constraint. Default value is `:none`.

### Exemple
```julia
julia > a, b = Variable("a", Complex), Variable("b", Complex)
julia > cstr = abs2(a^2) + abs2(3*b) << 25
julia > cstr.precond = :sqrt
```
"""
mutable struct Constraint
  p::Polynomial
  lb::Number
  ub::Number
  precond::Symbol

  Constraint(p::AbstractPolynomial, lb::Number, ub::Number) = new(p, lb, ub, :none)
end

#############################
## Constraint constructors
#############################
function <<(p::T, bound::Number) where T<: AbstractPolynomial
  return Constraint(p, -Inf-im*Inf, bound)
end

function <<(bound::Number, p::T) where T<: AbstractPolynomial
  return Constraint(p, bound, Inf+im*Inf)
end

function >>(p::T, bound::Number) where T<:AbstractPolynomial
  return bound << p
end

function >>(bound::Number, p::T) where T<:AbstractPolynomial
  return p << bound
end

function <<(bound::Number, cstr::Constraint)
  if real(bound) > real(cstr.ub) || imag(bound) > imag(cstr.ub)
    warn("<<(bnd, cstr): Creating a constraint with lower bound ", bound, " and upper bound ", cstr.ub)
  end
  return Constraint(cstr.p, bound, cstr.ub)
end
>>(bound::Number, cstr::Constraint) = cstr << bound

function <<(cstr::Constraint, bound::Number)
  if real(cstr.lb) > real(bound) || imag(cstr.lb) > imag(bound)
    warn("<<(cstr, bnd): Creating a constraint with lower bound ", cstr.lb, " and upper bound ", bound)
  end
  return Constraint(cstr.p, cstr.lb, bound)
end
>>(cstr::Constraint, bound::Number) = bound << cstr


function ==(p::T, x::Number) where T<:AbstractPolynomial
  return Constraint(p, x, x)
end


#############################
## Print
#############################
function Base.print(io::IO, cstr::Constraint)
  if cstr.lb == cstr.ub
    print(io, cstr.p, " = ", cstr.ub)
  else
    if cstr.lb != (-Inf-im*Inf)
      print(io, cstr.lb, " < ")
    end
    print(io, cstr.p)
    if cstr.ub !=(Inf+im*Inf)
      print(io, " < ", cstr.ub)
    end
  end
end