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
    C_N2::Float64
    T::Float64
end

struct AnodeMPLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    C_N2::Float64
    T::Float64
end

struct AnodeCLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    C_N2::Float64
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
    C_N2::Float64
    lambda::Float64
    T::Float64
    eta_c::Float64
end

struct CathodeMPLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    C_N2::Float64
    T::Float64
end

struct CathodeGDLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    C_N2::Float64
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

"""H2/O2/N2 species derivative contribution for porous layers and CLs."""
struct MEAGasSpeciesDerivative{nb_gdl, nb_mpl}
    agdl_C_H2::NTuple{nb_gdl, Float64}
    ampl_C_H2::NTuple{nb_mpl, Float64}
    acl_C_H2::Float64
    agdl_C_N2::NTuple{nb_gdl, Float64}
    ampl_C_N2::NTuple{nb_mpl, Float64}
    acl_C_N2::Float64
    ccl_C_O2::Float64
    cmpl_C_O2::NTuple{nb_mpl, Float64}
    cgdl_C_O2::NTuple{nb_gdl, Float64}
    ccl_C_N2::Float64
    cmpl_C_N2::NTuple{nb_mpl, Float64}
    cgdl_C_N2::NTuple{nb_gdl, Float64}
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




# ────────────────────────────────────────────────────────────────────────────────
# Assembly functions
# ────────────────────────────────────────────────────────────────────────────────

"""Assemble one complete MEA derivative from per-physics derivative contributions."""
function assemble_mea_derivative_1D(dw::MEADissolvedWaterDerivative,
                                    lw::MEALiquidWaterDerivative{NB_GDL, NB_MPL},
                                    vw::MEAVaporDerivative{NB_GDL, NB_MPL},
                                    sd::MEAGasSpeciesDerivative{NB_GDL, NB_MPL},
                                    vd::MEAVoltageDerivative,
                                    td::MEATemperatureDerivative{NB_GDL, NB_MPL}) where {NB_GDL, NB_MPL}

    # Gas-channel derivatives are still filled by gas-channel equations later.
    agc = AnodeGCDerivative(NaN, NaN, NaN, NaN, NaN)
    cgc = CathodeGCDerivative(NaN, NaN, NaN, NaN, NaN)

    agdl = ntuple(NB_GDL) do j
        AnodeGDLDerivative(vw.agdl_C_v[j], lw.agdl_s[j], sd.agdl_C_H2[j], sd.agdl_C_N2[j], td.agdl_T[j])
    end
    ampl = ntuple(NB_MPL) do j
        AnodeMPLDerivative(vw.ampl_C_v[j], lw.ampl_s[j], sd.ampl_C_H2[j], sd.ampl_C_N2[j], td.ampl_T[j])
    end
    acl = AnodeCLDerivative(vw.acl_C_v, lw.acl_s, sd.acl_C_H2, sd.acl_C_N2, dw.acl_lambda, td.acl_T)
    mem = MembraneDerivative(dw.mem_lambda, td.mem_T)
    ccl = CathodeCLDerivative(vw.ccl_C_v, lw.ccl_s, sd.ccl_C_O2, sd.ccl_C_N2, dw.ccl_lambda, td.ccl_T, vd.ccl_eta_c)
    cmpl = ntuple(NB_MPL) do j
        CathodeMPLDerivative(vw.cmpl_C_v[j], lw.cmpl_s[j], sd.cmpl_C_O2[j], sd.cmpl_C_N2[j], td.cmpl_T[j])
    end
    cgdl = ntuple(NB_GDL) do j
        CathodeGDLDerivative(vw.cgdl_C_v[j], lw.cgdl_s[j], sd.cgdl_C_O2[j], sd.cgdl_C_N2[j], td.cgdl_T[j])
    end

    return CellDerivative1D{NB_GDL, NB_MPL}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end

"""Assemble complete MEA derivatives by injecting GC contributions into existing MEA derivatives."""
function assemble_gc_derivative_1D(mea::AbstractVector{<:CellDerivative1D{NB_GDL, NB_MPL}},
                                   gas::GCGasDerivative{NB_GC},
                                   liq::GCLiquidWaterDerivative{NB_GC},
                                   temp::GCTemperatureDerivative{NB_GC}) where {NB_GDL, NB_MPL, NB_GC}
    return [begin
                agc = AnodeGCDerivative(gas.agc_C_v[i], liq.agc_s[i], gas.agc_C_H2[i], gas.agc_C_N2[i], temp.agc_T[i])
                cgc = CathodeGCDerivative(gas.cgc_C_v[i], liq.cgc_s[i], gas.cgc_C_O2[i], gas.cgc_C_N2[i], temp.cgc_T[i])
                node = mea[i]
                CellDerivative1D{NB_GDL, NB_MPL}(agc, node.agdl, node.ampl, node.acl, node.mem,
                                                    node.ccl, node.cmpl, node.cgdl, cgc)
            end for i in 1:NB_GC]
end
