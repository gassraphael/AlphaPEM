# cell_intermediates.jl
#
# Typed intermediate values used during the ODE right-hand side evaluation.
# Field sets are layer-specific to avoid carrying unused quantities.

abstract type AbstractCellIntermediate end

# Gas-channel intermediates need pressure/composition/transport properties.
struct AnodeGCIntermediates <: AbstractCellIntermediate
    P::Float64
    C_tot::Float64
    rho::Float64
    mu::Float64
end

struct CathodeGCIntermediates <: AbstractCellIntermediate
    P::Float64
    C_tot::Float64
    rho::Float64
    mu::Float64
end

# Thermal cache used by temperature dynamics inside the MEA.
struct MEAThermalIntermediates{nb_gdl, nb_mpl} <: AbstractCellIntermediate
    agdl::NTuple{nb_gdl, Float64}
    ampl::NTuple{nb_mpl, Float64}
    acl::Float64
    mem::Float64
    ccl::Float64
    cmpl::NTuple{nb_mpl, Float64}
    cgdl::NTuple{nb_gdl, Float64}
end

# Per-GC-column intermediates.
struct MEAIntermediates1D{nb_gdl, nb_mpl} <: AbstractCellIntermediate
    agc::AnodeGCIntermediates
    cgc::CathodeGCIntermediates
    rho_Cp0::MEAThermalIntermediates{nb_gdl, nb_mpl}
    i_n::Float64
    T_acl_mem_ccl::Float64
    v_re::Union{Nothing, Float64}
    k_purge::Union{Nothing, Float64}
end

# Manifold intermediates are represented even if some branches are temporarily disabled.
struct ManifoldNodeIntermediates <: AbstractCellIntermediate
    P::Float64
    Phi::Float64
    M::Union{Nothing, Float64}
    rho::Union{Nothing, Float64}
    mu::Union{Nothing, Float64}
end

struct ManifoldIntermediates{nb_man} <: AbstractCellIntermediate
    asm::NTuple{nb_man, ManifoldNodeIntermediates}
    aem::NTuple{nb_man, ManifoldNodeIntermediates}
    csm::NTuple{nb_man, ManifoldNodeIntermediates}
    cem::NTuple{nb_man, ManifoldNodeIntermediates}
end

struct AuxiliaryIntermediates <: AbstractCellIntermediate
    v_re::Union{Nothing, Float64}
    k_purge::Union{Nothing, Float64}
    Abp_a::Union{Nothing, Float64}
    Abp_c::Union{Nothing, Float64}
end

struct FuelCellIntermediatesP2D{nb_gdl, nb_mpl, nb_gc, MI, AI} <: AbstractCellIntermediate
    # MI is the manifold-intermediates container type (e.g. ManifoldIntermediates{nb_man}).
    # AI is the auxiliary-intermediates container type (e.g. AuxiliaryIntermediates).
    mea::NTuple{nb_gc, MEAIntermediates1D{nb_gdl, nb_mpl}}
    manifold::MI
    auxiliary::AI
end

