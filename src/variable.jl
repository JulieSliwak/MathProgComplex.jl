export Degree, Variable, show

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
## Type detection utils
#############################
isreal(x::Variable) = x.kind <: Real
iscomplex(x::Variable) = x.kind <: Complex
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