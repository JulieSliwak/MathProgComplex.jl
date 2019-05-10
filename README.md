# MathProgComplex.jl

Dev:
[![Build Status](https://travis-ci.org/JulieSliwak/MathProgComplex.jl.svg?branch=master)](https://travis-ci.org/JulieSliwak/MathProgComplex.jl)
[![codecov](https://codecov.io/gh/JulieSliwak/MathProgComplex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JulieSliwak/MathProgComplex.jl)

The `MathProgComplex` module is a tool for polynomial optimization problems with complex variables. These problems consist in optimizing a generic complex multivariate polynomial function, subject to some complex polynomial equality and inequality constraints.
The `MathProgComplex` module enables:

- the manipulation of multivariate polynomials with complex numbers to construct polynomial optimization problems with complex variables (POP-C).
- the evaluation of polynomials, for example the objective and the constraints of a (POP-C) from points
- the resolution of a (POP-C) via a JuMP model
- the export of a (POP-C) to be solved using another language

## Setting Julia for custom modules

The modules need to be accessible from julia's path to be loaded with the `using` command.
This can be done by running `push!(LOAD_PATH, "location/of/modules")`, which will update the path for the current session.

The path can be updated at every start of a Julia session by adding the command to the `.juliarc.jl` file, which should be located (or created) at the location given by `homedir()`.

## Structures

The `MathProgComplex` environment provides a structure and methods for working with **complex polynomial optimization problems** subject to **polynomial constraints**. It is based on a polynomial environment that allows to work on polynomial objects with natural operations (+, -, \*, conj, |.|^2).

### Polynomials

- The base type is `Variable`: it is a structure with a name (a string) and a type (Complex, Real or Bool).

```julia
using MathProgComplex
a = Variable("a", Complex)
b = Variable("b", Real)
c = Variable("c", Bool)
```

From `Variable` type, `Exponent` and `Polynomial` can be constructed by calling the respective constructors or with algebraic operations (+, -, \*, conj, |.|).

- An `Exponent` is a product of `Variables`.

```julia
expo1 = a*b
expo2 = Base.power_by_squaring(conj(a),3)*Base.power_by_squaring(b,5)
expo3 = abs2(a) # =a*conj(a)
```

> **_NOTE:_** `x^p` has been replaced by `Base.power_by_squaring(x,p)` in Julia 1.0.3


- A `Polynomial` is a sum of `Exponents` times complex numbers.

```julia
p = 3*expo1 + (4+2im)*expo2 +2im*expo3
```

- The `Point` type holds the variables at which polynomials can be evaluated.

#### Implemented methods

    - isconst, isone
    - evaluate
    - abs2, conj
    - is_homogeneous: tests if p(exp(iϕ)X) = p(X) ∀X∈R^n, Φ∈R (phaseinvariant  equation)
    - cplx2real: convert the provided object to a tuple of real and imaginary part, expressed with real and imaginary part variables.

```julia
using MathProgComplex

a = Variable("a", Complex)
b = Variable("b", Real)
p = a*conj(a) + b + 2

print(p)
# 2 + b + conj(a) * a

pt = Point([a, b], [1+2im, 1+im])
print(pt)
# a 1 + 2im
# b 1

pt = pt - Point([a], [1])
print(pt)
# a 0 + 2im
# b 1

val = evaluate(p, pt) # 7.0 + 0.0im

p_real, p_imag = cplx2real(p)
pt_r = cplx2real(pt)
val_real = evaluate(p_real, pt_r) # 7.0
val_imag = evaluate(p_imag, pt_r) # 0
```

### Polynomial optimization problems

- A `Constraint` structure holds a `Polynomial` and complex upper and lower bounds.

- A `Problem` is made up of:
  - the collection of `Variable`s of the problem (updated internally),
  - a `Polynomial` objective,
  - several named `constraint`s.

#### Implemented methods

- get_objective, set_objective!
- get_variables, get_variabletype, has_variable, add_variable!
- get_constraint, get_constraints, has_constraint, constraint_type, add_constraint!, rm_constraint!
- get_slacks, get_minslack
- pb_cplx2real: converts the problem variables, objective and constraints to real expressions function of real and imaginary part of the original problem variables.

Here is a full example, more can be found in the `examples` and `test` directories.

```julia
using MathProgComplex

a = Variable("a", Complex)
b = Variable("b", Real)
p_obj = abs2(a) + abs2(b) + 2
p_cstr1 = 3*a + b + 2
p_cstr2 = abs2(b) + 5*a*b + 2

pb = Problem()

set_objective!(pb, p_obj)

add_constraint!(pb, "Cstr 1", p_cstr1 << (3+5im)) # 2 + (3.0)*a + b < 3 + 5im

add_constraint!(pb, "Cstr 2", (2-im) << p_cstr2 << (3+7im)) # 2 - 1im < 2 + (5.0)*a * b + b^2 < 3 + 7im

add_constraint!(pb, "Cstr 3", p_cstr2 == 0) # 2 + (5.0)*a * b + b^2 = 0


print(pb)
# ▶ variables: a b
# ▶ objective: 2 + conj(a) * a + b^2
# ▶ constraints:
#  →     Cstr 1: 2 + (3.0)*a + b < 3 + 5im
#  →     Cstr 2: 2 - 1im < 2 + (5.0)*a * b + b^2 < 3 + 7im
#  →     Cstr 3: 2 + (5.0)*a * b + b^2 = 0

pt_sol = Point([a, b], [1, 1+2im])
print("slack by constraint at given point:\n", get_slacks(pb, pt_sol))
# slack by constraint at given point:
# Cstr 1 -3.0 + 5.0im
# Cstr 2 -5.0 + 1.0im
# Cstr 3 -8.0 - 0.0im
```

### Resolution

The polynomial optimization problems can be converted into JuMP models or be exported into formatted text files to be used in another language.

#### Using JuMP

 ```julia
m, JuMPvar = get_JuMP_cartesian_model(pb, solver)
solve(m)
 ```

#### Using file export

The `.dat` text format used in this module allows storing polynomial optimization problems in inary, Real or Complex variables, along with a scalar value for each variable (a point).

```julia
export_to_dat(pb, amplexportpath, point)
run_knitro(amplexportpath, amplscriptpath)
```
