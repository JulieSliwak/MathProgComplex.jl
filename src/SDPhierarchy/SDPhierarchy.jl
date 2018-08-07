# module SDPhierarchy

# using MathProgComplex
# using DataStructures
using OPFInstances

import JuMP, Mosek, MathProgBase, SCS, MathProgBase, CSDP

export RelaxationContext, Moment, MomentMatrix, SDPDual, SDPPrimal
export SDP_Instance, SDP_Block, SDP_CtrObjName, SDP_Problem


const DictType = SortedDict

###############################################################################
## Relaxation context, symmetries and cliques
###############################################################################
mutable struct RelaxationContext
    ismultiordered::Bool
    issparse::Bool
    symmetries::Set{DataType}       # ::SortedSet{DataType}
    hierarchykind::Symbol           # :Complex or :Real
    renamevars::Bool                # Replace variables with by shorter named ones
    di::Dict{String, Int}
    ki::Dict{String, Int}
    cstrtypes::Dict{String, Symbol}
    binvar_constraints::SortedDict{String, Constraint}      # equilibrium constraints x(1-x)=0 for modeling binary variables
    relaxparams::OrderedDict{Symbol, Any}

    RelaxationContext() = new(false,
                              false,
                              Set{DataType}(),
                              :Real,
                              false,
                              Dict{String, Int}(),
                              Dict{String, Int}(),
                              Dict{String, Symbol}(),
                              SortedDict{String, Constraint}(),
                              get_defaultparams())
end

include(joinpath("core", "build_relaxationcontext.jl"))

include(joinpath("core", "build_decomposition.jl"))

abstract type AbstractSymmetry end
type PhaseInvariance <: AbstractSymmetry end
include("symmetries.jl")

###############################################################################
## Moment Problem
###############################################################################
struct Moment
    conj_part::Exponent
    expl_part::Exponent
    clique::String
end

include(joinpath("base_types", "moment.jl"))

"""
    MomentMatrix{T}(mm, vars, order, matrixkind)

Store a moment or localizing matrix of size `order`, corresponding to the `vars` variables in the `mm` dictionnary.
**Note** that the matrix is indexed by a tuple of exponents, *the first of which contains only conjugated variables*, et second only real ones.
"""
mutable struct MomentMatrix{T}
    mm::DictType{Tuple{Exponent, Exponent}, DictType{Moment, T}}    # (row index, col index) -> lin. comb. of moments
    vars::Set{Variable}                                             # The set of variables from which the MomentMatrix was built
    order::Int                                                      # The order at which the moment matrix was built (d-ki for localizing constraints)
    matrixkind::Symbol                                              # Either :SDP, :SDPC or :Null
end

include(joinpath("base_types", "momentmatrix.jl"))

"""
    momentrel = SDPDual(obj, cstrs, moment_overlap)

Store a Moment Relaxation problem.
"""
struct SDPDual{T}
    objective::DictType{Moment, T}                                  # A linear comb. of moments, to be maximized
    constraints::DictType{Tuple{String, String}, MomentMatrix{T}}   # A set of moment matrices, either SDP or Null. A constraint (`key[1]`) can be split on several cliques (`key[2]`)
    moments_overlap::DictType{Exponent, Set{String}}                # A set of clique per exponent, describing coupling constraints
end

include(joinpath("core", "build_momentrelaxation.jl"))



###############################################################################
## SOS Problem
###############################################################################
const CtrName = Moment

mutable struct SDPPrimal{T}
    block_to_vartype::DictType{String, Symbol}                          # Either :SDP, :Sym, :SDPC, :SymC
    blocks::DictType{Tuple{CtrName, String, Exponent, Exponent}, T}     # (constraintname, block_name, γ, δ) -> coeff
    linsym::DictType{Tuple{CtrName, String, Exponent}, T}               # (constraintname, block_name, var) -> coeff
    lin::DictType{Tuple{CtrName, Exponent}, T}                          # (constraintname, var) -> coeff
    cst::DictType{CtrName, T}                                           #  constraintname -> coeff
end

include(joinpath("core", "build_SOSrelaxation.jl"))
include(joinpath("io", "export_SDPPrimal.jl"))


###############################################################################
## Solver structures
###############################################################################
type SDP_Instance
  VAR_TYPES
  BLOCKS
  LINEAR
  CONST
end


type SDP_Block
  id::Int64
  name::String
  var_to_id::SortedDict{String, Int64}

  SDP_Block(id::Int64, name::String) = new(id, name, SortedDict{String, Int64}())
end

const SDP_CtrObjName = Tuple{String, String, String}

"""
    SDP_Problem

Description of a SDP problem in the primal form, with string to integer maps coefficients matrices, scalar variables and contraint keys.

All SDP problems to be solved should be converted to this structure, for which the Mosek solver can be readily used, and can be easily extended to other solvers.

            max               ∑ A_0i[k,l] × Zi[k,l] + ∑ b_0[k] × x[k] + c_0
            s.t.    lb_j  <=  ∑ A_ji[k,l] × Zi[k,l] + ∑ b_j[k] × x[k] + c_j  <=  ub_j

**Notes**:
- Only the lower triangular part of coefficient matrices is stored. Hence a slice of the initial matrix is stroed, **no diagonal or non-diagonal coefficient is scaled**.
"""
type SDP_Problem
  # SDP vars
  name_to_sdpblock::SortedDict{String, SDP_Block}                                   # SDP variable name -> SDP variable description
  id_to_sdpblock::SortedDict{Int64, SDP_Block}                                      # SDP variable id   -> SDP variable description

  # Scalar variables
  scalvar_to_id::Dict{String, Int64}                                                # Scalar variable name -> scalar variable id

  # Objective / constraints
  obj_keys::SortedSet{SDP_CtrObjName}                                               # Set of SDP_CtrObjName corresponding to the objective
  name_to_ctr::SortedDict{SDP_CtrObjName, Tuple{Int64, String, Float64, Float64}}   # SDP_CtrObjName constraint -> constraint id j, type, lower and upper bounds (lb_j, ub_j)
  id_to_ctr::SortedDict{Int64, SDP_CtrObjName}                                      # SDP_CtrObjName constraint id -> SDP_CtrObjName name

  matrices::SortedDict{Tuple{SDP_CtrObjName, String, String, String}, Float64}      # matrix coefficients           (j, i, k, l) -> A_ji[k,l]
  linear::SortedDict{Tuple{SDP_CtrObjName, String}, Float64}                        # scalar variables coeffs       (j, k) -> b_j[k]
  cst_ctr::SortedDict{SDP_CtrObjName, Float64}                                      # constant coeff                j -> c_j

  SDP_Problem() = new(SortedDict{String, SDP_Block}(),
                      SortedDict{Int64, SDP_Block}(),
                      Dict{String, Int64}(),
                      SortedSet{SDP_CtrObjName}(),
                      SortedDict{SDP_CtrObjName, Tuple{Int64, String, Float64, Float64}}(),
                      SortedDict{Int64, SDP_CtrObjName}(),
                      SortedDict{Tuple{SDP_CtrObjName, String, String, String}, Float64}(),
                      SortedDict{Tuple{SDP_CtrObjName, String}, Float64}(),
                      SortedDict{SDP_CtrObjName, Float64}()
                      )
end

include(joinpath("solvers", "Mosek.jl"))
include(joinpath("solvers", "JuMP.jl"))


include(joinpath("SDP_Instance", "common.jl"))
include(joinpath("SDP_Instance", "build_from_sdpfile.jl"))
include(joinpath("SDP_Instance", "build_from_SDPPrimal.jl"))
include(joinpath("SDP_Instance", "build_from_SDPDual.jl"))
# include(joinpath("SDP_Instance", "fileexport.jl"))


###############################################################################
## Unsorted
###############################################################################


include(joinpath("core", "run_hierarchy.jl"))

include("example_problems.jl")
include("utils.jl")
include(joinpath("io", "momentsos_io.jl"))
include(joinpath("io", "print.jl"))

# end