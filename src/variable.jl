export Degree, Variable, show, iscomplex, isreal, isbool

"""
    Degree(explvar::Int, conjvar::Int)

Define a mathematical degree, which is used to define a variable exponent to
`explvar` and its conjugate exponent to `conjvar`.

### Attributes
- `explvar`: Int
- `conjvar`: Int
"""
mutable struct Degree
    explvar::Int
    conjvar::Int
end

#############################
## Equality, hash
#############################
==(d1::Degree, d2::Degree) = (d1.explvar==d2.explvar) && (d1.conjvar==d2.conjvar)
!=(d1::Degree, d2::Degree) = !(d1 == d2)
hash(deg::Degree, h::UInt) = hash(deg.explvar, hash(deg.conjvar, h))



"""
Variable(varname::String, Complex)

Define a mathematical variable `varname` of a certain mathematical kind,
`Complex` here.

### Attributes
- `name`: String
- `kind`: a type among `Complex`, `Real` and `Bool`
"""
mutable struct Variable <: AbstractPolynomial
    name::String
    kind::Type

    function Variable(name, kind)
        validtypes = Set([Complex, Real, Bool])
        if true ∉ map(x-> kind <: x, validtypes) || typeof(name) ∉ Set([String, SubString{String}])
            error("Variable() : attempting to define a variable $name of type $kind, supported types are {Complex, Real, Bool}")
        end
        return new(String(name), kind)
    end
end

#############################
## Equality, hash
#############################
==(x::Variable, y::Variable) = (x.name==y.name) && (x.kind==y.kind)
!=(x::Variable, y::Variable) = !(x == y)
hash(var::Variable, h::UInt) = hash(var.name, hash(var.kind, h))


#############################
## Type detection utils
#############################
iscomplex(x::Variable) = x.kind <: Complex
isreal(x::Variable) = x.kind <: Real
isbool(x::Variable) = x.kind <: Bool

#############################
## Show
#############################
function Base.show(io::IO, x::Variable)
    print(io, x.name)
end

function Base.show(io::IO, d::Degree)
    print(io, "(", d.explvar, ",", d.conjvar, ")")
end