# MathProgComplex.jl

The `MathProgComplex` module enables the construction of polynomial optimization problems in complex variables for power systems. It uses the `PolynomialOptim` environment. The environment can adapt to several input formats (Matpower, GOC, IIDM) if the user defines the proper functions (a read function and elements of the power network along with their parameters and constraints). Once the optimization problem is built, the problem can be exported in real numbers and in text files (to be treated by AMPL for example). It is also possible to convert the problem in a JuMP model in real numbers (cartesian or polar form).

<!-- ## Setting Julia for custom modules
The modules need to be accessible from julia's path to be loaded with the `using` command.
This can be done by running `push!(LOAD_PATH, "location/of/modules")`, which will update the path for the current session.

The path can be updated at every start of a Julia session by adding the command to the `.juliarc.jl` file, which should be located (or created) at the location given by `homedir()`. -->

## PolynomialOptim

The `PolynomialOptim` environment provides a structure and methods for working with **complex polynomial optimization problems** subject to **polynomial constraints**. It is based on a polynomial environment that allows to work on polynomial objects with natural operations (+, -, \*, conj, |.|).

The base type is `Variable`, from which `Exponents`, `Monomial`, and `Polynomial` can be constructed by calling the respective constructors or with algebraic operations (+, -, \*, conj, |.|). The `Point` type holds the variables at which polynomials can be evaluated.

Available functions on `Polynomial`, `Mononomial`, `Variable` types:

- isconst, isone
- evaluate
- abs2, conj
- is_homogeneous: tests if p(exp(iϕ)X) = p(X) ∀X∈R^n, Φ∈R (phase invariant equation)
- cplx2real: convert the provided object to a tuple of real and imaginary part, expressed with real and imaginary part variables.

```julia
include(joinpath("src_PolynomialOptim", "PolynomialOptim.jl"))

a = Variable("a", Complex)
b = Variable("b", Real)
p = a*conj(a) + b + 2

print(p)
# (2.0) + b + conj(a) * a

pt = Point([a, b], [1+2im, 1+im])
print(pt)
# a 1 + 2im
# b 1

pt = pt - Point([a], [1])
print(pt)
# a 0 + 2im
# b 1

val = evaluate(p, pt) # 7 + 0im

p_real, p_imag = cplx2real(p)
pt_r = cplx2real(pt)
val_real = evaluate(p_real, pt_r) # 7 + 0im
val_imag = evaluate(p_imag, pt_r) # 0
```

### Polynomial problems

A `Constraint` structure holding a `Polynomial` and complex upper and lower bounds is defined to build the `Problem` type made of:

- several `Variables`
- a `Polynomial` objective,
- several constraints names along with their corresponding `Constraint`

### Implemented methods

- get_objective, set_objective!
- get_variables, get_variabletype, has_variable, add_variable!
- get_constraint, get_constraints, has_constraint, constraint_type, add_constraint!, rm_constraint!
- get_slacks, get_minslack
- pb_cplx2real: converts the problem variables, objective and constraints to real expressions function of real and imaginary part of the original problem variables.

```julia
include(joinpath("src_PolynomialOptim", "PolynomialOptim.jl"))

a = Variable("a", Complex)
b = Variable("b", Real)
p_obj = abs2(a) + abs2(b) + 2
p_cstr1 = 3a + b + 2
p_cstr2 = abs2(b) + 5a*b + 2

pb = Problem()

add_variable!(pb, a); add_variable!(pb, b)

set_objective!(pb, p_obj)

add_constraint!(pb, "Cstr 1", p_cstr1 << 3+5im) # -∞ ≤ (2.0) + b^2 + (5.0) * a * b ≤ 2 + 6im
add_constraint!(pb, "Cstr 2", 2-im << p_cstr2 << 3+7im) # -∞ ≤ (2.0) + b^2 + (5.0) * a * b ≤ 2 + 6im
add_constraint!(pb, "Cstr 3", p_cstr2 == 0) # -∞ ≤ (2.0) + b^2 + (5.0) * a * b ≤ 2 + 6im

print(pb)
# ▶ variables: a b
# ▶ objective: (2.0) + conj(a) * a + b^2
# ▶ constraints:
#     Cstr 1: -∞ < (2.0) + (3.0) * a + b < 3 + 5im
#     Cstr 2: 2 - 1im < (2.0) + (5.0) * a * b + b^2 < 3 + 7im
#     Cstr 3: (2.0) + (5.0) * a * b + b^2 = 0

pt_sol = Point([a, b], [1, 1+2im])
print("slack by constraint at given point:\n", get_slacks(pb, pt_sol))
# slack by constraint at given point:
# Cstr 1 -3.0 + 0.0im
# Cstr 2 -5 + 1im
# Cstr 3 -8 + 0im
```
