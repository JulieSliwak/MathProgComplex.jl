# MathProgComplex.jl Documentation

A julia package for polynomial optimization in complex numbers.

## Package features

The package is comprised of two submodules:

- PolynomialOptim offers:
  - Variables, exponents, polynomials structures and natural algebra operations,
  - Problem and constraint structures for polynomial optimization under polynomial constraints (POP), along with functions for evaluating a candidate point wrt the current problem,
  - export of polynomial problems to JuMP models, for subsequent solve with any non linear solver,
  - export to a formatted text file, so as to be able to work outside of the Julia JuMP environment
- SDPhierarchy offers enables to build and solve SDP relaxations of POP, the so-called Lasserre hierarchy :
  - high level workflow for building and solving SDP relaxations of input POP with given parameters
  - structures for primal and dual sdps
  - interface to Mosek using Mosek.jl C API, and interface to JuMP models
  - file export fomat of SDP problems

## Package outline

[package outline]