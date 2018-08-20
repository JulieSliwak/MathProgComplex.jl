# Guide

This page aims to present the main features of the `SDPhierarchy` submodule, and both high-level and low-level workflows. Details on how to set the hierarchy high-level parameters is available at the page [Relaxation options](@ref), a mathematic description of each problem along with important implementation details is available at [Implementation details](@ref).

## Purpose

This module aims to implement a hierarchy of *Semi-Definite Positive* relaxations for Polynomials Optimization Problems ($POP$), modeled with the `PolynomialOptim` sub-module of MathProgComplex. Given an input $POP$, one can choose to build a more or less tight SDP relaxation having certain properties, leading to more or less large and numerically tractable problems. This hierarchy was first laid out by J.B. Lasserre for real polynomial optimization problems. It was extended to complex polynomial optimization problems by C. Josz and D. Molzahn in order to leverage the underlying structure present in complex polynomial problems, lost when converted to euivalent real polynomial problems by expliciting real and imaginary parts. Such complex polynomial optimization problems, $POP-\mathbb C$ arise for example when working on AC models of current transmission networks ($ACOPF$) - for which current at each node of the network is modeled as a complex variable.

The SDP hierarchy applied to $POP-\mathbb R$, the $POP-\mathbb C$ converted to real variables, yields exact relaxations with orders of either 1 or 2 empirically for $ACPOF$. However, when dealing with problems having hundreds or more variables, the order 2 SDP relaxation is too large to be handled and solved by current state-of-the-art SDP solvers. Hence the will to use the structure present in complex problems, exploit symmetries, specify relaxation orders constraint-wise (and not problem-wise) and decompose variables.

### References

- Global optimization with polynomials and the problem of moments, *Lasserre, Jean B*, SIAM Journal on optimization, 2001;
- Moment/sum-of-squares hierarchy for complex polynomial optimization, *Josz, Cédric and Molzahn, Daniel K*, arXiv preprint arXiv:1508.02068, 2015;
- Multi-ordered Lasserre hierarchy for large scale polynomial optimization in real and complex variables, *Josz, Cedric and Molzahn, Daniel K*, arXiv preprint arXiv:1709.04376, 2017.

## Features of the hierarchy

The following list details the main features of the implemented SDP hierarchy:

- *real or complex hierarchy*: given a $POP-\mathbb C$, one can build a complex SDP relaxation, generally more tractable than its real counterpart. Ideally, this problem would be solved by complex SDP solvers, however we don't know of any existing complex SDP solvers mature enough at this date. Therefore the complex SDP is converted to a real SDP and then solved, which is still generally preferable to switching to real numbers at the POP stage. Of course, the real hierarchy can also be constructed from the $POP-\mathbb R$.
- *multi-ordered*: one $POP$ can lead to many more or less tights relaxations, depending on a global relaxation order $d$ [^1]. Actually, one can choose the relaxation orders at a constraint level. This feature can yield intermediate sized SDPs assuming a good decomposition is used. More details are available at section [Maximal cliques description](@ref).
- *dense or sparse*: given a $POP$ with a relaxation order $d$, the resulting SDP problem will have a matrix variable with as many rows and columns as exponents of global degree up to $d$ that can be formed from the variables of the $POP$. Therefore the size of SDO relaxations grow exponentially with the global order $d$ (or maximal constraint order $d_i$). However one can build decomposed versions of this problem, allowing to split this large SDP matrix into smaller ones, along with coupling constraints. Finding the best decomposition given an input $POP$ and constraint wise relaxation orders is a difficult open problem, not tackled by this module.
- *symmetries*: it is known that if the input $POP$ objective and constraints polynomials show a certain symmetry, all non-null moments of the SDP relaxation solutions will also respect this symmetry, which allows to set some coefficients of the SDP variables to 0 from the begining. One symmetry arising on current transmission network is *phase invariance*: $\forall z\in \mathbb C, \theta\in [0, 2\pi[,\quad p(z) = p(e^{i \theta}z)$.
- *binary variables*: a variable $b\in\{0,1\}$ appearing in an optimization problem can be treated as a continuous variable $x^b$ with the added polynomial constraint $x^b (1-x^b)=0$, which fits well the framework of the POP hierarchy.

[^1]: The global degree of an exponent $\bar{z}^\alpha z^\beta = \prod_i \bar{z}^{\alpha_i} z^{\beta_i}$ is $d = |\alpha| + |\beta| = \sum_i \alpha_i + \beta_i$.

## Hierarchy workflow

The hierarchy workflow is as follows.

Expected input is :

- the $POP$: $\min\left\{ f(x) : g_i(x) \ge 0,  i=1, ..., m\right\}$,
- a `RelaxationContext` object containing all relaxation choices,
- a `max_cliques` object, describing the decomposition chosen for the relaxation. For a dense relaxation, it can be built with `get_maxcliques(relax_ctx, problem)`.

Then several steps are taken:

- Parameters defining the constraints of the SDP problem are computed. Specifically, they give the orders $d_i$, $d_i - k_i$ and relevant variables for the moment matrix $\mathcal M_{d_i^{Cl_i}}(y_i)$ and localizing matrices $\mathcal M_{d_i-k_i}(g_i y_i)$. See [Building the moment relaxation](@ref) for more details.
- The moment relaxation of order $(d_1, \ldots, d_m)$, stored as a `SDPDual` object is computed,
- Its lagrangian dual, the sum of squares (SOS) relaxation, stored as a `SDPPrimal` can be computed, possibly stored to a formatted file format,
- either one of the relaxation problems can be converted into a `SDP_Instance`, which holds the SDP problem as a primal problem, with integer indexation for solvers,
- The `SDP_Instance` can be passed to the Mosek solver with the C API, or converted into a JuMP problem, that can be solved with any JuMP(compatible SDP solver (Mosek, CSDP, SCS, ...).

![global_workflow](../images/hierarchy_workflow.png)

## Simple "high level" workflow

All parameters are hold by the `RelaxationContext` object, the previously described workflow is executed upon calling `run_hierarchy`. Options are detailed at section [Relaxation options](@ref)

```julia
## Build polynomial problem
x1 = Variable("x", Real)
x2 = Variable("y", Real)
problem = Problem()
set_objective!(problem, -1.0*x1)
add_constraint!(problem, "ineq", (x1^2+x2^2) << 4)
θ1 = π/3
add_constraint!(problem, "eq_rot1", (cos(θ1)*x1+sin(θ1)*x2) == 0)

## Set parameters
relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                    d = 1,
                                    params = Dict(:opt_outlev=>1,
                                                  :opt_pbsolved=>:SOSRelaxation,
                                                  :opt_solver=>:MosekCAPI))

## Build the order 1 SOS relaxation and pass it to Mosek via the C API.
primobj, dualobj = run_hierarchy(problem, relax_ctx);
```

## Simple "lower level" example

Here is the low level equivalent of the previous code example.

```julia
using MathProgComplex

## Build polynomial problem
x1 = Variable("x", Real)
x2 = Variable("y", Real)
problem = Problem()
set_objective!(problem, -1.0*x1)
add_constraint!(problem, "ineq", (x1^2+x2^2) << 4)
add_constraint!(problem, "eq_rot1", (cos(π/3)*x1+sin(π/3)*x2) == 0)

relax_ctx = set_relaxation(problem; hierarchykind=:Real,
                                    d = 1,
                                    params = Dict(:opt_outlev=>1))


# Build sparsity pattern, chordal extension, maximal cliques.
# Simple dense problem : one clique with all variables.
max_cliques = get_maxcliques(relax_ctx, problem)

# Compute moment and localizing matrices parameters: order et variables
momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

# Build the moment relaxation problem
momentrel = build_momentrelaxation(relax_ctx, problem, momentmat_param,
                                        localizingmat_param, max_cliques)

# Convert to a primal SDP problem
sosrel = build_SOSrelaxation(relax_ctx, momentrel)

# sdp_instance_moment::SDP_Problem = build_SDP_Instance_from_SDPDual(momentrel)
sdp_instance_sos::SDP_Problem = build_SDP_Instance_from_SDPPrimal(sosrel)

primal = SortedDict{Tuple{String,String,String}, Float64}()
dual = SortedDict{Tuple{String, String, String}, Float64}()

primobj, dualobj = solve_JuMP(sdp_instance_sos, :CSDPSolver, primal, dual;
                                              sol_info = relax_ctx.relaxparams,
                                              optsense = :Max)
```