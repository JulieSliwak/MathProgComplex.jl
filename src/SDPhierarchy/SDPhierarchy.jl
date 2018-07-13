# module SDPhierarchy

# using MathProgComplex
# using DataStructures
using OPFInstances

export RelaxationContext, Moment, MomentMatrix, SDPDual, SDPPrimal
export SDP_Instance, SDP_Block, SDP_Moment, SDP_Problem


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
    relaxparams::OrderedDict{Symbol, Any}

    RelaxationContext() = new(false,
                              false,
                              Set{DataType}(),
                              :Real,
                              false,
                              Dict{String, Int}(),
                              Dict{String, Int}(),
                              Dict{String, Symbol}(),
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
    mm::DictType{Tuple{Exponent, Exponent}, DictType{Moment, T}}
    vars::Set{Variable}
    order::Int
    matrixkind::Symbol            # Either :SDP or :Sym
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
mutable struct SDPPrimal{T}
    block_to_vartype::DictType{String, Symbol}                       # Either :SDP, :Sym, :SDPc, :SymC
    blocks::DictType{Tuple{Moment, String, Exponent, Exponent}, T}   # ((α, β), block_name, γ, δ) -> coeff
    linsym::DictType{Tuple{Moment, String, Exponent}, T}             # ((α, β), block_name, var) -> coeff
    lin::DictType{Tuple{Moment, Exponent}, T}                        # ((α, β), var) -> coeff
    cst::DictType{Moment, T}                                         # (α, β) -> coeff
end

include(joinpath("core", "build_SOSrelaxation.jl"))
include(joinpath("io", "export_SDPPrimal.jl"))
include("SDPPrimal_cplx2real.jl")



###############################################################################
## Mosek Structures
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

const SDP_Moment = Tuple{String, String, String}

type SDP_Problem
  # SDP vars
  name_to_sdpblock::SortedDict{String, SDP_Block}
  id_to_sdpblock::SortedDict{Int64, SDP_Block}

  # Scalar variables
  scalvar_to_id::Dict{String, Int64}

  # Objective / constraints
  obj_keys::SortedSet{SDP_Moment}
  name_to_ctr::SortedDict{SDP_Moment, Tuple{Int64, String, Float64, Float64}} # Id, type et bornes des contraintes
  id_to_ctr::SortedDict{Int64, SDP_Moment}

  matrices::SortedDict{Tuple{SDP_Moment, String, String, String}, Float64} # Matrices SDP du corps des contraintes / objectif
  linear::SortedDict{Tuple{SDP_Moment, String}, Float64} # Matrice portant les parties linéaires des contraintes
  cst_ctr::SortedDict{SDP_Moment, Float64} # Constante du corps des contraintes

  SDP_Problem() = new(SortedDict{String, SDP_Block}(),
                      SortedDict{Int64, SDP_Block}(),
                      Dict{String, Int64}(),
                      SortedSet{SDP_Moment}(),
                      SortedDict{SDP_Moment, Tuple{Int64, String, Float64, Float64}}(),
                      SortedDict{Int64, SDP_Moment}(),
                      SortedDict{Tuple{SDP_Moment, String, String, String}, Float64}(),
                      SortedDict{Tuple{SDP_Moment, String}, Float64}(),
                      SortedDict{SDP_Moment, Float64}()
                      )
end

include(joinpath("solvers", "run_mosek.jl"))

include(joinpath("SDP_Instance", "common.jl"))
include(joinpath("SDP_Instance", "build_from_sdpfile.jl"))
include(joinpath("SDP_Instance", "build_from_SDPPrimal.jl"))



###############################################################################
## Unsorted
###############################################################################


include(joinpath("core", "run_hierarchy.jl"))

include("example_problems.jl")
include("utils.jl")
include(joinpath("io", "momentsos_io.jl"))
include(joinpath("io", "print.jl"))

# end