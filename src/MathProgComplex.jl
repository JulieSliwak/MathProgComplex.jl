isdefined(Base, :__precompile__) && __precompile__()

module MathProgComplex

using DataStructures, JuMP

# see https://invenia.github.io/Memento.jl/latest/faq/pkg-usage.html#Using-Memento-in-Julia-packages?-1
using Memento
using Compat: @__MODULE__
const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

# message hierarchy: debug < info < warn < error
# do not print efficiency-related warnings
setlevel!(LOGGER, "info")

# print efficiency-related warnings
# setlevel!(LOGGER, "debug")


import Base: ==, !=, <<, >>, isless, isconst, isreal, isnull, isequal
import Base: +, -, *, /, ^, conj, conj!, abs2, norm, real, imag
import Base: show, print, convert, copy, hash, merge
import Base: start, next, done, length, setindex!, getindex, haskey, keys, values, deepcopy

export AbstractPolynomial

abstract type  AbstractPolynomial end

include("utils_internal.jl")

include("variable.jl")
include("exponent.jl")
include("polynomial.jl")

include("constraint.jl")
include("problem.jl")
include("point.jl")

include("cplx2real.jl")
include("evaluate.jl")
include("iterators.jl")
include("utils_Poly.jl")
include("problem_accessors.jl")

## algebra
include(joinpath("algebra", "add_algebra.jl"))
include(joinpath("algebra", "mult_algebra.jl"))
include(joinpath("algebra", "order.jl"))
include(joinpath("algebra", "unaries.jl"))

## export dat
include(joinpath("export_dat", "utils_ampl.jl"))
include(joinpath("export_dat", "utils_dat_export.jl"))
include(joinpath("export_dat", "utils_dat_import.jl"))

## export JuMP
include(joinpath("export_JuMP", "utils_jump.jl"))


# SDPhierarchy function
include(joinpath("SDPhierarchy", "SDPhierarchy.jl"))

end
