# Implementation details

This sections aims at explaining the choices of structures and functions. We recall the global workflow:

![global_workflow](../images/hierarchy_workflow.png)

## Input data

- `problem::Problem`
- `max_cliques::SortedDict{String, Set{Variable}}`
- `relax_ctx::RelaxationContext`

### The $POP$

The $POP-\mathbb R$ and $POP-\mathbb C$:

```math
\begin{aligned}
\min_{x\in\mathbb R^n} & f(x) &\\
\text{s.t.} & gi(x) \ge 0 & i=1,\ldots, m\\
\end{aligned}
\qquad
\begin{aligned}
\min_{z\in\mathbb C^n} & f(z) &\\
\text{s.t.} & gi(z) \ge 0 & i=1,\ldots, m\\
\end{aligned}
```

### Maximal cliques description

From the input $POP$, a sparsity pattern can be constructed. It is the graph having one node per variable, and one edge for each pair of variables appearing in the same exponent. Given the relaxation order $d_i$ and the constraint polynomial degree $k_i$ of each constraint, if $d_i > k_i$, then edges have to be added to represent all the pair of variables appearing in the considered *constraint*.

Given this extended sparsity pattern, any chordal extension will provide a decomposition for the considered relaxation. One builds a chordal extension of a given graph by adding edges so that all cycles have length 4 or more have at least a chord (an edge between two non-consecutive nodes of the cycle).

The second expected input is the maximal clique decomposition of this chordal graph. It allows to decompose one large SDP matrix into several smaller ones, tied with coupling constraints. For dense decompositions, only one clique containing all variables should be provided. It can by constructed by calling `get_maxcliques`:

```julia
max_cliques = get_maxcliques(relax_ctx, problem)
```

See *Moment/sum-of-squares hierarchy for complex polynomial optimization*, p.27 for an example.

### `RelaxationContext`

The `RelaxationContext` serves as a general structure for holding the user choices, some pre-computed information such as the relaxation, polynomial degree and constraint type of each constraint, or added constraints for binary variables of the input $POP$. Besides it also holds the `relaxparams` ordered dictionary, that serves to hold all timings and measures made while running the hierarchy, printed in the final csv. It is constructed by calling `set_relaxation`.

## Building the moment relaxation

The goal is to build the $(d_1, \ldots, d_m)$ moment relaxation, with the maximal clique decomposition $C_1, \ldots, C_p$.

At this point some parameters have to be computed :

- the minimal union of cliques $I_i$ necessary to have all variables appearing in $g_i$, and therefore form the localizing matrix $\mathcal M_{d_i-k_i}(g_i y, I_i)$,
- the minimum order $d_i^{cl}$ possible for the moment matrix associated with clique $C_i$, given the order of the constraints that have variables in clique $i$^.

See *Moment/sum-of-squares hierarchy for complex polynomial optimization*, p.18 for details.

Then, the moment relaxation can be derived:

```math
\begin{aligned}
\inf_{y\in \mathcal H} & \sum_{\alpha, \beta} f_{\alpha, \beta} y_{\alpha, \beta} & \\
\text{s.t.}            & y_{0,0} = 1 & \\
                       & \mathcal M_{d_l^{cl}}(y, C_l) \succeq 0, & l = 1, \ldots, p, \\
                       & \mathcal M_{d_i-k_i}(g_i y, I_i) \succeq 0, & l = 1, \ldots, m, \\
\end{aligned}
```

In this context of decomposed problems, a moment $y_{\alpha, \beta, l}$ is uniquely defined by an exponent $\bar{z}^\alpha z^\beta$, and the clique $C_l$ from which its variables are taken from. Indeed, as soon as decomposition is used, there are moments that can be expressed in several cliques. The simplest of which are the ones representing variables that belong to several cliques.

The matrices $\mathcal M_{d_l^{cl}}(y, C_l)$ have elements at all indexes $(\bar{z}^\gamma, z^\delta)$, where $|\gamma|, |\delta| \le d_l^{cl}$ and $z$ has vars in clique $C_l$. The element at index $(\bar{z}^\gamma, z^\delta)$ of this matrix is the moment $y_{\gamma, \delta, l}$ associated with clique $l$.

The matrices $\mathcal M_{d_i-k_i}(g_i y, I_i)$ have elements at all indexes $(\bar{z}^\gamma, z^\delta)$, where $|\gamma|, |\delta| \le d_i-k_i$ and $z$ has vars in variables set $I_i$. The element at index $(\bar{z}^\gamma, z^\delta)$ of this matrix is the linear combination of moments $\sum_{\alpha, \beta} g_{\alpha, \beta} y_{\alpha+\gamma, \beta+\delta, l}$ associated with one of the cliques making up $I_i$ and comprising all variables from exponent $\bar{z}^{\alpha+\delta} z^{\beta+\gamma}$. Currently, one clique is arbitrarily chosen.

### moment relaxation parameters

The following parameters are obtained by calling `build_sparsity(relax_ctx, problem, max_cliques)`:

- the integer $d_i^{cl}$ per clique $i$ to build the moment matrix $\mathcal M_{d_i^{cl}}(y, C_i)$
- the integer $d_i - k_i$ and the set of cliques $L$ required to build $I_i = \bigcup_{l\in L} C_l$.

```julia
momentmat_param, localizingmat_param = build_sparsity(relax_ctx, problem, max_cliques)

momentmat_param::OrderedDict{String, Int}
localizingmat_param::OrderedDict{String, Tuple{Set{String}, Int}}
```

### Moment relaxation

The moment relaxation is stored by a `SDPDual` object:

```julia
struct SDPDual{T}
    objective::DictType{Moment, T}                                  # A linear comb. of moments, to be maximized
    constraints::DictType{Tuple{String, String}, MomentMatrix{T}}   # A set of moment matrices, either SDP or Null. A constraint (`key[1]`) can be split on several cliques (`key[2]`)
    moments_overlap::DictType{Exponent, Set{String}}                # A set of clique per exponent, describing coupling constraints
end
```

It is assembled as follows:

1. First the objective of the $POP$ $\sum_{\alpha, \beta} f_{\alpha, \beta} \bar{z}^\alpha z^\beta$ is converted to a linear combination of moment $\sum_{\alpha, \beta} f_{\alpha, \beta} y_{\alpha, \beta, clique_i}$. Each moment has to be associated with a clique that fully supports it.

2. Each clique $C_l$ result in a `MomentMatrix`, representing $\mathcal M_{d_l^{cl}}(y, C_l)$. Here all moments are attached to the clique $C_l$.

3. Each one-sided inequality of the input $POP$ is then put on the form $\sum_{\alpha, \beta} g_{i, \alpha, \beta} \bar{z}^\alpha z^\beta \ge 0$, from which a `MomentMatrix` representing $\mathcal M_{d_i-k_i}(g_i y, I_i)$ is built. Double sided inequalities are split into two one sided inequalities.
    Equality constraint $g_i(z)$ can either be treated as two constraints $g_i(z) \le 0$ and $g_i(z) \ge 0$, or as imposing the localizing matrix $\mathcal M_{d_i-k_i}(g_i y, I_i)$ to be null.

4. Finally, if decomposition is used, the list of exponents that can be expressed in at least two cliques is computed, in order to lay out the right coupling constraint.

## Building the SOS relaxation

The SOS relaxation is stored by a `SDPPrimal` object:

```julia
mutable struct SDPPrimal{T}
    block_to_vartype::DictType{String, Symbol}                          # Either :SDP, :Sym, :SDPC, :SymC
    blocks::DictType{Tuple{CtrName, String, Exponent, Exponent}, T}     # (constraintname, block_name, γ, δ) -> coeff
    linsym::DictType{Tuple{CtrName, String, Exponent}, T}               # (constraintname, block_name, var) -> coeff
    lin::DictType{Tuple{CtrName, Exponent}, T}                          # (constraintname, var) -> coeff
    cst::DictType{CtrName, T}                                           #  constraintname -> coeff
end
```

It is built from a `SDPDual` instance, as its lagrangian dual.

## Storing and solving SDPs

In order to separate the process of building SDP problems and solving them, a structure `SDP_Problem` has been defined. It carries a problem indexed with string keys to keep readability, and mappings from string to integers for passing to solvers. The problem is expected to be in the following form:

```math
\begin{aligned}
& \min_{Z_i\succeq 0, x_i\in \mathbb{R}}
& \sum_i A_{0i} \cdot Z_i + \sum_i b_{0i} x_i + c_0 \\
& \text{s.t.}
& \sum_i A_{ji} \cdot Z_i + \sum_i b_{ji} x_i + c_j =0\\
&&Z_i \succeq 0, \; i=1,\ldots, n_{sdp}&\\
&&x_i\in\mathbb{R}, \; i=1,\ldots,n_{scal}
\end{aligned}
```

Both `SDPPrimal` and `SDPDual` can be converted into `SDP_Problem` instances. This instances can be saved to text files with a homemade file format `*.sdp`, and imported back. Finally, these problems can be solved by using the C API of Mosek, defined in the package `Mosek.jl`, or by converting them to JuMP models, that can then be solved by any SDP solvers. The high level workflow allows to work with Mosek, CSDP and SCS. Experiments have shown CSDP to be reliable, as opposed to SCS.

!!! note

    There is a bug in JuMP, concerning maximization problems having a non null constant term in their objective, see [JuMP issue 1390](https://github.com/JuliaOpt/JuMP.jl/issues/1390). This renders the package tests using the open source CSDP incorrect.

## Possible improvements

- symmetries: all well with SOS relaxation, as less constraints, but for moment problem, are solvers able to leverage the fact that coefficients from the SDP matrix variable are null ? The phase invariance symmetry applied to real $ACOPF$ relaxations should yield moment matrices that are quite sparse, but not in an orderly fashion. For complex SDP relaxation, the moment matrix variable(s) would have a nice block diagonal structure. If solvers can leverage this sparsity *on the variables*, and not the problem definition as they usually do, important improvements in memory usage (sparse storage for some of the vectors) and CPU time (cholesky factorizations of sparse block diagonal matrices, and fast sparse matrix matrix products) are to be expected.

- data structure: dictionaries are extensively used throughout this module, either the `Base` `dict` or the `DataStructures.SortedDict` and `DataStructures.OrderedDict`. The question of using the best structure at each step of the workflow has not been clearly answered. There is much to gain in clarifying this matter, at least to make the implementation more robust, and to gain performance. Right now the dictionary type for many structures is an internal variable of the SDPhierarchy module `DictType`.