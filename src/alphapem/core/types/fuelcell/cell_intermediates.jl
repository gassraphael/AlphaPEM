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

"""Reusable vector workspace shared by `calculate_flow_1D_GC_manifold_int_values!`
and `velocity_profiles_from_inlet_flows!`.

Both functions compute the same GC thermodynamic quantities (pressures, mole
fractions, viscosities), so a single pre-allocated workspace covers both call
sites and prevents any heap allocation.  Fields specific to one function are
documented inline.
"""
mutable struct GCManifoldWorkspace
    # ── Shared: GC thermodynamics (filled by both velocity and manifold functions) ──
    P_agc       :: Vector{Float64}   # Anode GC ideal-gas pressure (Pa)
    P_cgc       :: Vector{Float64}   # Cathode GC ideal-gas pressure (Pa)
    Phi_agc     :: Vector{Float64}   # Anode GC relative humidity (–)
    Phi_cgc     :: Vector{Float64}   # Cathode GC relative humidity (–)
    M_agc       :: Vector{Float64}   # Anode GC mixture molar mass (kg·mol⁻¹)
    M_cgc       :: Vector{Float64}   # Cathode GC mixture molar mass (kg·mol⁻¹)
    rho_agc     :: Vector{Float64}   # Anode GC mixture density (kg·m⁻³)
    rho_cgc     :: Vector{Float64}   # Cathode GC mixture density (kg·m⁻³)
    y_H2_agc    :: Vector{Float64}   # Dry H₂ mole fraction in anode GC
    y_O2_cgc    :: Vector{Float64}   # Dry O₂ mole fraction in cathode GC
    x_H2O_v_agc :: Vector{Float64}   # H₂O vapour mole fraction in anode GC
    x_H2_agc    :: Vector{Float64}   # H₂ mole fraction in anode GC
    x_N2_agc    :: Vector{Float64}   # N₂ mole fraction in anode GC
    x_H2O_v_cgc :: Vector{Float64}   # H₂O vapour mole fraction in cathode GC
    x_O2_cgc    :: Vector{Float64}   # O₂ mole fraction in cathode GC
    x_N2_cgc    :: Vector{Float64}   # N₂ mole fraction in cathode GC
    mu_gaz_agc  :: Vector{Float64}   # Anode GC dynamic viscosity (Pa·s)
    mu_gaz_cgc  :: Vector{Float64}   # Cathode GC dynamic viscosity (Pa·s)
    # ── Velocity-specific: interface fluxes and flow/pressure profiles ───────────
    C_tot_agdl     :: Vector{Float64}   # Total molar concentration at AGC/AGDL interface
    C_tot_cgdl     :: Vector{Float64}   # Total molar concentration at CGDL/CGC interface
    J_tot_agc_agdl :: Vector{Float64}   # Net molar flux AGC→AGDL (mol·m⁻²·s⁻¹)
    J_tot_cgdl_cgc :: Vector{Float64}   # Net molar flux CGDL→CGC (mol·m⁻²·s⁻¹)
    J_a         :: Vector{Float64}   # Anode channel molar flow along the GC (mol·m⁻²·s⁻¹)
    J_c         :: Vector{Float64}   # Cathode channel molar flow along the GC (mol·m⁻²·s⁻¹)
    P_a_chan    :: Vector{Float64}   # Anode channel pressure profile with viscous drop (Pa)
    P_c_chan    :: Vector{Float64}   # Cathode channel pressure profile with viscous drop (Pa)
    v_a         :: Vector{Float64}   # Anode gas velocity (m·s⁻¹)
    v_c         :: Vector{Float64}   # Cathode gas velocity (m·s⁻¹)
    v_a_nominal :: Vector{Float64}   # Anode velocity reordered for counter-flow
end

function GCManifoldWorkspace(nb_gc::Int)
    z = () -> Vector{Float64}(undef, nb_gc)
    return GCManifoldWorkspace(
        z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(),  # shared GC thermodynamics (12)
        z(), z(), z(), z(), z(), z(),                                 # mole fractions (x_H2O_v_agc … x_N2_cgc)
        z(), z(), z(), z(), z(), z(), z(), z(), z(), z(), z(),       # velocity-specific (11)
    )
end

