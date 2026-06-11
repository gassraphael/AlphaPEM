# cell_derivative.jl
#
# Typed representation of time derivatives used in the ODE right-hand side.
# These structures mirror the state topology from cell_state.jl.

abstract type AbstractCellDerivative end

# Anode side
struct AnodeGCDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    C_N2::Float64
    T::Float64
end

struct AnodeGDLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    T::Float64
end

struct AnodeMPLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    T::Float64
end

struct AnodeCLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    lambda::Float64
    T::Float64
end

# Electrolyte
struct MembraneDerivative <: AbstractCellDerivative
    lambda::Float64
    T::Float64
end

# Cathode side
struct CathodeCLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    lambda::Float64
    T::Float64
    eta_c::Float64
end

struct CathodeMPLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    T::Float64
end

struct CathodeGDLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    T::Float64
end

struct CathodeGCDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    C_N2::Float64
    T::Float64
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold derivatives
# ────────────────────────────────────────────────────────────────────────────────

"""Time derivatives for pressure and relative humidity at a manifold node."""
struct ManifoldDerivative <: AbstractCellDerivative
    P::Float64      # dP/dt
    Phi::Float64    # dPhi/dt
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D cell-column derivative (MEA + AGC/CGC)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete 1D derivative container for one gas-channel column."""
struct CellDerivative1D{nb_gdl, nb_mpl}
    agc::AnodeGCDerivative
    agdl::NTuple{nb_gdl, AnodeGDLDerivative}
    ampl::NTuple{nb_mpl, AnodeMPLDerivative}
    acl::AnodeCLDerivative
    mem::MembraneDerivative
    ccl::CathodeCLDerivative
    cmpl::NTuple{nb_mpl, CathodeMPLDerivative}
    cgdl::NTuple{nb_gdl, CathodeGDLDerivative}
    cgc::CathodeGCDerivative
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D MEA-core derivative contributions (per-physics, before final assembly)
# ────────────────────────────────────────────────────────────────────────────────

"""Dissolved-water contribution (lambda only) for ACL, membrane and CCL."""
struct MEADissolvedWaterDerivative
    acl_lambda::Float64
    mem_lambda::Float64
    ccl_lambda::Float64
end

"""Liquid-water contribution (saturation derivatives) across porous layers."""
struct MEALiquidWaterDerivative{nb_gdl, nb_mpl}
    agdl_s::NTuple{nb_gdl, Float64}
    ampl_s::NTuple{nb_mpl, Float64}
    acl_s::Float64
    ccl_s::Float64
    cmpl_s::NTuple{nb_mpl, Float64}
    cgdl_s::NTuple{nb_gdl, Float64}
end

"""Water-vapour contribution (C_v derivatives) across porous layers and CLs."""
struct MEAVaporDerivative{nb_gdl, nb_mpl}
    agdl_C_v::NTuple{nb_gdl, Float64}
    ampl_C_v::NTuple{nb_mpl, Float64}
    acl_C_v::Float64
    ccl_C_v::Float64
    cmpl_C_v::NTuple{nb_mpl, Float64}
    cgdl_C_v::NTuple{nb_gdl, Float64}
end

"""H2/O2 species derivative contribution for porous layers and CLs."""
struct MEAH2O2SpeciesDerivative{nb_gdl, nb_mpl}
    agdl_C_H2::NTuple{nb_gdl, Float64}
    ampl_C_H2::NTuple{nb_mpl, Float64}
    acl_C_H2::Float64
    ccl_C_O2::Float64
    cmpl_C_O2::NTuple{nb_mpl, Float64}
    cgdl_C_O2::NTuple{nb_gdl, Float64}
end

"""Voltage contribution (eta_c derivative in the CCL)."""
struct MEAVoltageDerivative
    ccl_eta_c::Float64
end

"""Temperature contribution across porous layers, CLs and membrane."""
struct MEATemperatureDerivative{nb_gdl, nb_mpl}
    agdl_T::NTuple{nb_gdl, Float64}
    ampl_T::NTuple{nb_mpl, Float64}
    acl_T::Float64
    mem_T::Float64
    ccl_T::Float64
    cmpl_T::NTuple{nb_mpl, Float64}
    cgdl_T::NTuple{nb_gdl, Float64}
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D gas-channel derivative contributions (per-physics, before final assembly)
# ────────────────────────────────────────────────────────────────────────────────

"""Gas-species contribution across anode/cathode gas channels for all GC nodes."""
struct GCGasDerivative{nb_gc}
    agc_C_v::NTuple{nb_gc, Float64}
    agc_C_H2::NTuple{nb_gc, Float64}
    agc_C_N2::NTuple{nb_gc, Float64}
    cgc_C_v::NTuple{nb_gc, Float64}
    cgc_C_O2::NTuple{nb_gc, Float64}
    cgc_C_N2::NTuple{nb_gc, Float64}
end

"""Liquid-water saturation contribution across anode/cathode gas channels for all GC nodes."""
struct GCLiquidWaterDerivative{nb_gc}
    agc_s::NTuple{nb_gc, Float64}
    cgc_s::NTuple{nb_gc, Float64}
end

"""Temperature contribution across anode/cathode gas channels for all GC nodes."""
struct GCTemperatureDerivative{nb_gc}
    agc_T::NTuple{nb_gc, Float64}
    cgc_T::NTuple{nb_gc, Float64}
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold line derivatives
# ────────────────────────────────────────────────────────────────────────────────

"""Complete derivative container for one manifold line."""
struct ManifoldLineDerivative{nb_nodes}
    nodes::NTuple{nb_nodes, ManifoldDerivative}
end

# ────────────────────────────────────────────────────────────────────────────────
# P2D fuel-cell derivatives (cell columns = MEA + AGC/CGC)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete P2D derivative container for the cell-column stack.
"""
struct FuelCellDerivativeP2D{nb_gdl, nb_mpl, nb_gc}
    nodes::NTuple{nb_gc, CellDerivative1D{nb_gdl, nb_mpl}}
end


