# SOS hierarchy

So far, the ``src_SOShierarchy`` code is able to trnasform a ``Problem`` object into a description of a primal SDP problem corresponding to a specified SOS hierarchy relaxation.

## Current workflow

```julia
    # Construction of the initial problem
    problem = buildPOP_WB2_expl()

    # Set relaxation parameters
    relax_ctx = set_relaxation(problem, hierarchykind=:Complex, d = 2, leveragesymmetries = true)

    # Build sparsity pattern, compute maximal cliques
    max_cliques = get_maxcliques(relax_ctx, problem)

    # Compute moment matrices parameters: order et variables
    moments_params = build_sparsity(relax_ctx, problem, max_cliques)

    # Compute partial moment hierarchy
    mmtrel_pb = MomentRelaxationPb(relax_ctx, problem, moments_params, max_cliques)

    # Convert to a primal SDP problem
    SDP_body, SDP_rhs = build_SDP(relax_ctx, mmtrel_pb)
```

## Relaxation parameters

- ``:Real`` or ``:Complex`` : according to the input problem kind, build the corresponding hierarchy.
- Constraint specific relaxation order (`ismultiordered`)
- Global or constraint specific relaxation order specification (`d` or `di` parameter)
- Can leverage symmetries in the initial problem to reduce the number of moments to be optimized over. Currently, one symmetry is implemented : $f(ze^{i\theta}) = f(z) \quad \forall z\in\mathbb C, \theta\in[0,2\pi[$. All moments associated with monomials $\bar z^\alpha z^\beta$ such that $|\alpha| \neq |\beta|$ are fixed to 0.

## Relevant features

## To Be Done

**13.04**:

- Enforce triang storage everywhere
- no symmetries for this current example ? (real)