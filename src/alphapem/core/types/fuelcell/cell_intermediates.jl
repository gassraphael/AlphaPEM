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
struct CellIntermediates1D{nb_gdl, nb_mpl} <: AbstractCellIntermediate
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
    mea::NTuple{nb_gc, CellIntermediates1D{nb_gdl, nb_mpl}}
    manifold::MI
    auxiliary::AI
end


# ────────────────────────────────────────────────────────────────────────────────
# Reusable vector workspaces for intermediate calculations
# ────────────────────────────────────────────────────────────────────────────────

"""Reusable vector workspace for `calculate_heat_int_values!`.
Holds only the effective thermal conductivity vectors that are otherwise
reallocated every call.
"""
mutable struct MEAHeatIntWorkspace
    k_th_eff_agdl_agdl::Vector{Float64}
    k_th_eff_ampl_ampl::Vector{Float64}
    k_th_eff_cmpl_cmpl::Vector{Float64}
    k_th_eff_cgdl_cgdl::Vector{Float64}
end

function MEAHeatIntWorkspace(nb_gdl::Int, nb_mpl::Int)
    return MEAHeatIntWorkspace(
        Vector{Float64}(undef, max(nb_gdl - 1, 0)),
        Vector{Float64}(undef, max(nb_mpl - 1, 0)),
        Vector{Float64}(undef, max(nb_mpl - 1, 0)),
        Vector{Float64}(undef, max(nb_gdl - 1, 0))
    )
end

"""Reusable vector workspace for `calculate_flows_1D_MEA_int_values!`.
Holds the capillary and diffusion coefficient vectors that are otherwise
reallocated every call to `flows_1D_MEA_int_values`.
"""
mutable struct MEAFlowsIntWorkspace
    D_cap_agdl_agdl::Vector{Float64}
    D_cap_ampl_ampl::Vector{Float64}
    D_cap_cmpl_cmpl::Vector{Float64}
    D_cap_cgdl_cgdl::Vector{Float64}
    Da_eff_agdl_agdl::Vector{Float64}
    Da_eff_ampl_ampl::Vector{Float64}
    Dc_eff_cmpl_cmpl::Vector{Float64}
    Dc_eff_cgdl_cgdl::Vector{Float64}
end

function MEAFlowsIntWorkspace(nb_gdl::Int, nb_mpl::Int)
    return MEAFlowsIntWorkspace(
        Vector{Float64}(undef, max(nb_gdl - 1, 0)),
        Vector{Float64}(undef, max(nb_mpl - 1, 0)),
        Vector{Float64}(undef, max(nb_mpl - 1, 0)),
        Vector{Float64}(undef, max(nb_gdl - 1, 0)),
        Vector{Float64}(undef, max(nb_gdl - 1, 0)),
        Vector{Float64}(undef, max(nb_mpl - 1, 0)),
        Vector{Float64}(undef, max(nb_mpl - 1, 0)),
        Vector{Float64}(undef, max(nb_gdl - 1, 0))
    )
end

"""Reusable vector workspace for `calculate_flow_1D_GC_manifold_int_values!`.
Holds the per-GC-node vectors that are otherwise reallocated every call.
"""
mutable struct GCManifoldWorkspace
    P_agc::Vector{Float64}
    P_cgc::Vector{Float64}
    Phi_agc::Vector{Float64}
    Phi_cgc::Vector{Float64}
    y_H2_agc::Vector{Float64}
    y_O2_cgc::Vector{Float64}
    M_agc::Vector{Float64}
    M_cgc::Vector{Float64}
    rho_agc::Vector{Float64}
    rho_cgc::Vector{Float64}
    mu_gaz_agc::Vector{Float64}
    mu_gaz_cgc::Vector{Float64}
end

function GCManifoldWorkspace(nb_gc::Int)
    return GCManifoldWorkspace(
        Vector{Float64}(undef, nb_gc), Vector{Float64}(undef, nb_gc),
        Vector{Float64}(undef, nb_gc), Vector{Float64}(undef, nb_gc),
        Vector{Float64}(undef, nb_gc), Vector{Float64}(undef, nb_gc),
        Vector{Float64}(undef, nb_gc), Vector{Float64}(undef, nb_gc),
        Vector{Float64}(undef, nb_gc), Vector{Float64}(undef, nb_gc),
        Vector{Float64}(undef, nb_gc), Vector{Float64}(undef, nb_gc)
    )
end

