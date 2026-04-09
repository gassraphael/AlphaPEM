# cell_flows.jl
#
# Typed representation of inter-layer fluxes and source terms produced during
# the ODE right-hand side evaluation.

# ────────────────────────────────────────────────────────────────────────────────
# Water management – scalar flux / source terms
# ────────────────────────────────────────────────────────────────────────────────

"""Sorption (absorption) source terms for water in the catalyst layers. Units: mol·m⁻³·s⁻¹
"""
struct MEASorptionSources # S_abs
    v_acl :: Float64   # Vapour absorption at the anode CL
    l_acl :: Float64   # Liquid absorption at the anode CL
    v_ccl :: Float64   # Vapour absorption at the cathode CL
    l_ccl :: Float64   # Liquid absorption at the cathode CL
end


"""Dissolved-water fluxes through the ionomer / membrane. Units: mol·m⁻²·s⁻¹
"""
struct MEADissolvedWaterFlux # J_lambda
    acl_mem :: Float64   # ACL → membrane flux
    mem_ccl :: Float64   # Membrane → CCL flux
end


"""Water production source terms at the catalyst layers (reaction + crossover). Units: mol·m⁻³·s⁻¹
"""
struct MEAWaterProductionSources # Sp
    acl :: Float64   # Water produced at the anode CL
    ccl :: Float64   # Water produced at the cathode CL
end


# ────────────────────────────────────────────────────────────────────────────────
# Gas species – scalar reaction / crossover source terms
# ────────────────────────────────────────────────────────────────────────────────

"""Gas reaction and crossover source terms (H₂ or O₂). Units: mol·m⁻³·s⁻¹
"""
struct MEAGasReactionSources
    reac :: Float64   # Electrochemical reaction consumption
    cros :: Float64   # Crossover loss
end


# ────────────────────────────────────────────────────────────────────────────────
# Thermal – scalar heat source terms
# ────────────────────────────────────────────────────────────────────────────────

"""Reaction heat sources at the catalyst layers. Units: W·m⁻³
"""
struct MEAReactionHeat # Q_r
    acl :: Float64
    ccl :: Float64
end


"""Sorption (absorption / desorption) heat sources at the catalyst layers. Units: W·m⁻³
"""
struct MEASorptionHeat # Q_sorp
    v_acl :: Float64
    l_acl :: Float64
    v_ccl :: Float64
    l_ccl :: Float64
end


"""Protonic (ionic) Joule-heating sources. Units: W·m⁻³
"""
struct MEAProtonHeat # Q_p
    mem :: Float64
    ccl :: Float64
end


# ────────────────────────────────────────────────────────────────────────────────────
# Through-plane inter-layer fluxes  (dimension-dependent)
#
# Type parameters used throughout:
#   NB_GDL : number of GDL nodes per side
#   NB_MPL : number of MPL nodes per side
#
# Inter-node interface counts are derived in place as `NB_GDL - 1` and `NB_MPL - 1`.
# ────────────────────────────────────────────────────────────────────────────────────

"""Liquid-water inter-layer fluxes through the MEA (anode + cathode porous layers). Units: kg·m⁻²·s⁻¹
"""
struct MEALiquidFluxes{NB_GDL, NB_MPL} # Jl
    agc_agdl  :: Float64
    agdl_agdl :: Vector{Float64}   # nb_gdl - 1 inter-node flows
    agdl_ampl :: Float64
    ampl_ampl :: Vector{Float64}   # nb_mpl - 1 inter-node flows
    ampl_acl  :: Float64
    ccl_cmpl  :: Float64
    cmpl_cmpl :: Vector{Float64}
    cmpl_cgdl :: Float64
    cgdl_cgdl :: Vector{Float64}
    cgdl_cgc  :: Float64

    # Constructor checks mesh-size invariants for inter-node arrays (NB_* - 1).
    # Inputs are normalized to Vector{Float64} for predictable downstream typing.
    function MEALiquidFluxes{NB_GDL, NB_MPL}(agc_agdl, agdl_agdl, agdl_ampl, ampl_ampl,
                                             ampl_acl, ccl_cmpl, cmpl_cmpl, cmpl_cgdl,
                                             cgdl_cgdl, cgdl_cgc) where {NB_GDL, NB_MPL}
        length(agdl_agdl) == NB_GDL - 1 || throw(ArgumentError("agdl_agdl length must be NB_GDL - 1."))
        length(ampl_ampl) == NB_MPL - 1 || throw(ArgumentError("ampl_ampl length must be NB_MPL - 1."))
        length(cmpl_cmpl) == NB_MPL - 1 || throw(ArgumentError("cmpl_cmpl length must be NB_MPL - 1."))
        length(cgdl_cgdl) == NB_GDL - 1 || throw(ArgumentError("cgdl_cgdl length must be NB_GDL - 1."))
        return new{NB_GDL, NB_MPL}(agc_agdl, collect(Float64, agdl_agdl), agdl_ampl,
                                   collect(Float64, ampl_ampl), ampl_acl, ccl_cmpl,
                                   collect(Float64, cmpl_cmpl), cmpl_cgdl,
                                   collect(Float64, cgdl_cgdl), cgdl_cgc)
    end
end

"""Water-vapour inter-layer fluxes through the MEA (anode + cathode porous layers). Units: mol·m⁻²·s⁻¹
"""
struct MEAVaporFluxes{NB_GDL, NB_MPL} # Jv
    agc_agdl  :: Float64
    agdl_agdl :: Vector{Float64}
    agdl_ampl :: Float64
    ampl_ampl :: Vector{Float64}
    ampl_acl  :: Float64
    ccl_cmpl  :: Float64
    cmpl_cmpl :: Vector{Float64}
    cmpl_cgdl :: Float64
    cgdl_cgdl :: Vector{Float64}
    cgdl_cgc  :: Float64

    # Constructor checks mesh-size invariants for inter-node arrays (NB_* - 1).
    # Inputs are normalized to Vector{Float64} for predictable downstream typing.
    function MEAVaporFluxes{NB_GDL, NB_MPL}(agc_agdl, agdl_agdl, agdl_ampl, ampl_ampl,
                                            ampl_acl, ccl_cmpl, cmpl_cmpl, cmpl_cgdl,
                                            cgdl_cgdl, cgdl_cgc) where {NB_GDL, NB_MPL}
        length(agdl_agdl) == NB_GDL - 1 || throw(ArgumentError("agdl_agdl length must be NB_GDL - 1."))
        length(ampl_ampl) == NB_MPL - 1 || throw(ArgumentError("ampl_ampl length must be NB_MPL - 1."))
        length(cmpl_cmpl) == NB_MPL - 1 || throw(ArgumentError("cmpl_cmpl length must be NB_MPL - 1."))
        length(cgdl_cgdl) == NB_GDL - 1 || throw(ArgumentError("cgdl_cgdl length must be NB_GDL - 1."))
        return new{NB_GDL, NB_MPL}(agc_agdl, collect(Float64, agdl_agdl), agdl_ampl,
                                   collect(Float64, ampl_ampl), ampl_acl, ccl_cmpl,
                                   collect(Float64, cmpl_cmpl), cmpl_cgdl,
                                   collect(Float64, cgdl_cgdl), cgdl_cgc)
    end
end

"""Hydrogen inter-layer fluxes (anode side only). Units: mol·m⁻²·s⁻¹
"""
struct MEAHydrogenFluxes{NB_GDL, NB_MPL} # J_H2
    agc_agdl  :: Float64
    agdl_agdl :: Vector{Float64}
    agdl_ampl :: Float64
    ampl_ampl :: Vector{Float64}
    ampl_acl  :: Float64

    # Constructor checks mesh-size invariants for inter-node arrays (NB_* - 1).
    # Inputs are normalized to Vector{Float64} for predictable downstream typing.
    function MEAHydrogenFluxes{NB_GDL, NB_MPL}(agc_agdl, agdl_agdl, agdl_ampl,
                                               ampl_ampl, ampl_acl) where {NB_GDL, NB_MPL}
        length(agdl_agdl) == NB_GDL - 1 || throw(ArgumentError("agdl_agdl length must be NB_GDL - 1."))
        length(ampl_ampl) == NB_MPL - 1 || throw(ArgumentError("ampl_ampl length must be NB_MPL - 1."))
        return new{NB_GDL, NB_MPL}(agc_agdl, collect(Float64, agdl_agdl), agdl_ampl,
                                   collect(Float64, ampl_ampl), ampl_acl)
    end
end

"""Oxygen inter-layer fluxes (cathode side only). Units: mol·m⁻²·s⁻¹
"""
struct MEAOxygenFluxes{NB_GDL, NB_MPL} # J_O2
    ccl_cmpl  :: Float64
    cmpl_cmpl :: Vector{Float64}
    cmpl_cgdl :: Float64
    cgdl_cgdl :: Vector{Float64}
    cgdl_cgc  :: Float64

    # Constructor checks mesh-size invariants for inter-node arrays (NB_* - 1).
    # Inputs are normalized to Vector{Float64} for predictable downstream typing.
    function MEAOxygenFluxes{NB_GDL, NB_MPL}(ccl_cmpl, cmpl_cmpl, cmpl_cgdl,
                                             cgdl_cgdl, cgdl_cgc) where {NB_GDL, NB_MPL}
        length(cmpl_cmpl) == NB_MPL - 1 || throw(ArgumentError("cmpl_cmpl length must be NB_MPL - 1."))
        length(cgdl_cgdl) == NB_GDL - 1 || throw(ArgumentError("cgdl_cgdl length must be NB_GDL - 1."))
        return new{NB_GDL, NB_MPL}(ccl_cmpl, collect(Float64, cmpl_cmpl), cmpl_cgdl,
                                   collect(Float64, cgdl_cgdl), cgdl_cgc)
    end
end


# ────────────────────────────────────────────────────────────────────────────────
# Per-layer scalar source terms  (dimension-dependent)
# ────────────────────────────────────────────────────────────────────────────────

"""Liquid-water phase-change source terms at each porous-layer node. Units: mol·m⁻³·s⁻¹
"""
struct MEALiquidSources{NB_GDL, NB_MPL} # Sl
    agdl :: NTuple{NB_GDL, Float64}
    ampl :: NTuple{NB_MPL, Float64}
    acl  :: Float64
    ccl  :: Float64
    cmpl :: NTuple{NB_MPL, Float64}
    cgdl :: NTuple{NB_GDL, Float64}
end

"""Water-vapour phase-change source terms at each porous-layer node. Units: mol·m⁻³·s⁻¹
"""
struct MEAVaporSources{NB_GDL, NB_MPL} # Sv
    agdl :: NTuple{NB_GDL, Float64}
    ampl :: NTuple{NB_MPL, Float64}
    acl  :: Float64
    ccl  :: Float64
    cmpl :: NTuple{NB_MPL, Float64}
    cgdl :: NTuple{NB_GDL, Float64}
end


# ────────────────────────────────────────────────────────────────────────────────
# Thermal inter-layer fluxes and per-layer heat sources  (dimension-dependent)
# ────────────────────────────────────────────────────────────────────────────────

"""Thermal (conductive) inter-layer fluxes through the full MEA stack.
Two extra fields (`acl_mem`, `mem_ccl`) compared to the mass-flux structs.
Units: W·m⁻²
"""
struct MEAThermalFluxes{NB_GDL, NB_MPL} # Jt
    agc_agdl  :: Float64
    agdl_agdl :: Vector{Float64}
    agdl_ampl :: Float64
    ampl_ampl :: Vector{Float64}
    ampl_acl  :: Float64
    acl_mem   :: Float64   # ACL → membrane interface
    mem_ccl   :: Float64   # Membrane → CCL interface
    ccl_cmpl  :: Float64
    cmpl_cmpl :: Vector{Float64}
    cmpl_cgdl :: Float64
    cgdl_cgdl :: Vector{Float64}
    cgdl_cgc  :: Float64

    # Constructor checks mesh-size invariants for inter-node arrays (NB_* - 1).
    # Inputs are normalized to Vector{Float64} for predictable downstream typing.
    function MEAThermalFluxes{NB_GDL, NB_MPL}(agc_agdl, agdl_agdl, agdl_ampl, ampl_ampl,
                                              ampl_acl, acl_mem, mem_ccl, ccl_cmpl,
                                              cmpl_cmpl, cmpl_cgdl, cgdl_cgdl,
                                              cgdl_cgc) where {NB_GDL, NB_MPL}
        length(agdl_agdl) == NB_GDL - 1 || throw(ArgumentError("agdl_agdl length must be NB_GDL - 1."))
        length(ampl_ampl) == NB_MPL - 1 || throw(ArgumentError("ampl_ampl length must be NB_MPL - 1."))
        length(cmpl_cmpl) == NB_MPL - 1 || throw(ArgumentError("cmpl_cmpl length must be NB_MPL - 1."))
        length(cgdl_cgdl) == NB_GDL - 1 || throw(ArgumentError("cgdl_cgdl length must be NB_GDL - 1."))
        return new{NB_GDL, NB_MPL}(agc_agdl, collect(Float64, agdl_agdl), agdl_ampl,
                                   collect(Float64, ampl_ampl), ampl_acl, acl_mem,
                                   mem_ccl, ccl_cmpl, collect(Float64, cmpl_cmpl),
                                   cmpl_cgdl, collect(Float64, cgdl_cgdl), cgdl_cgc)
    end
end

"""Liquefaction / evaporation heat sources at each porous-layer node. Units: W·m⁻³
"""
struct MEALiquidHeat{NB_GDL, NB_MPL} # Q_liq
    agdl :: NTuple{NB_GDL, Float64}
    ampl :: NTuple{NB_MPL, Float64}
    acl  :: Float64
    ccl  :: Float64
    cmpl :: NTuple{NB_MPL, Float64}
    cgdl :: NTuple{NB_GDL, Float64}
end

"""Electronic Joule-heating sources at each porous-layer node. Units: W·m⁻³
"""
struct MEAElectricHeat{NB_GDL, NB_MPL} # Q_e
    agdl :: NTuple{NB_GDL, Float64}
    ampl :: NTuple{NB_MPL, Float64}
    acl  :: Float64
    ccl  :: Float64
    cmpl :: NTuple{NB_MPL, Float64}
    cgdl :: NTuple{NB_GDL, Float64}
end


# ────────────────────────────────────────────────────────────────────────────────
# Top-level flow containers  (future return types of the calculation functions)
# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
# Along-channel (GC-level) flow structures
#
# These cover the fluxes that travel along the gas channel direction (x-direction)
# between GC nodes, as well as the GC inlet/outlet boundary fluxes.
# GC-to-MEA interface fluxes (e.g., agc_agdl) are stored per-node in MEAFlows1D.
#
# Type parameter:
#   NB_GC : number of gas-channel nodes (spatial discretisation along the channel).
# ────────────────────────────────────────────────────────────────────────────────

"""Along-channel water-vapour molar fluxes at the GC level. Units: mol·m⁻²·s⁻¹

Fields agc_agc and cgc_cgc are inter-node fluxes of length NB_GC - 1.
"""
struct GCVaporFlows{NB_GC}
    agc_in  :: Float64          # Inlet vapour flux at the anode GC entry
    agc_agc :: Vector{Float64}  # Along-channel inter-node vapour fluxes (length NB_GC - 1)
    agc_out :: Float64          # Outlet vapour flux at the anode GC exit
    cgc_in  :: Float64          # Inlet vapour flux at the cathode GC entry
    cgc_cgc :: Vector{Float64}  # Along-channel inter-node vapour fluxes (length NB_GC - 1)
    cgc_out :: Float64          # Outlet vapour flux at the cathode GC exit

    function GCVaporFlows{NB_GC}(agc_in, agc_agc, agc_out,
                                  cgc_in, cgc_cgc, cgc_out) where {NB_GC}
        length(agc_agc) == NB_GC - 1 || throw(ArgumentError("agc_agc length must be NB_GC - 1 = $(NB_GC - 1)."))
        length(cgc_cgc) == NB_GC - 1 || throw(ArgumentError("cgc_cgc length must be NB_GC - 1 = $(NB_GC - 1)."))
        return new{NB_GC}(Float64(agc_in), collect(Float64, agc_agc), Float64(agc_out),
                          Float64(cgc_in), collect(Float64, cgc_cgc), Float64(cgc_out))
    end
end

"""Along-channel liquid-water mass fluxes at the GC level. Units: kg·m⁻²·s⁻¹

Liquid water has no inlet in either gas channel (only exits).
Fields agc_agc and cgc_cgc are inter-node fluxes of length NB_GC - 1.
"""
struct GCLiquidFlows{NB_GC}
    agc_agc :: Vector{Float64}  # Along-channel inter-node liquid fluxes (length NB_GC - 1)
    agc_out :: Float64          # Outlet liquid flux at the anode GC exit
    cgc_cgc :: Vector{Float64}  # Along-channel inter-node liquid fluxes (length NB_GC - 1)
    cgc_out :: Float64          # Outlet liquid flux at the cathode GC exit

    function GCLiquidFlows{NB_GC}(agc_agc, agc_out, cgc_cgc, cgc_out) where {NB_GC}
        length(agc_agc) == NB_GC - 1 || throw(ArgumentError("agc_agc length must be NB_GC - 1 = $(NB_GC - 1)."))
        length(cgc_cgc) == NB_GC - 1 || throw(ArgumentError("cgc_cgc length must be NB_GC - 1 = $(NB_GC - 1)."))
        return new{NB_GC}(collect(Float64, agc_agc), Float64(agc_out),
                          collect(Float64, cgc_cgc), Float64(cgc_out))
    end
end

"""Along-channel hydrogen molar fluxes at the anode GC level. Units: mol·m⁻²·s⁻¹

Field agc_agc is the inter-node flux vector of length NB_GC - 1.
"""
struct GCHydrogenFlows{NB_GC}
    agc_in  :: Float64          # Inlet H₂ flux at the anode GC entry
    agc_agc :: Vector{Float64}  # Along-channel inter-node H₂ fluxes (length NB_GC - 1)
    agc_out :: Float64          # Outlet H₂ flux at the anode GC exit

    function GCHydrogenFlows{NB_GC}(agc_in, agc_agc, agc_out) where {NB_GC}
        length(agc_agc) == NB_GC - 1 || throw(ArgumentError("agc_agc length must be NB_GC - 1 = $(NB_GC - 1)."))
        return new{NB_GC}(Float64(agc_in), collect(Float64, agc_agc), Float64(agc_out))
    end
end

"""Along-channel oxygen molar fluxes at the cathode GC level. Units: mol·m⁻²·s⁻¹

Field cgc_cgc is the inter-node flux vector of length NB_GC - 1.
"""
struct GCOxygenFlows{NB_GC}
    cgc_in  :: Float64          # Inlet O₂ flux at the cathode GC entry
    cgc_cgc :: Vector{Float64}  # Along-channel inter-node O₂ fluxes (length NB_GC - 1)
    cgc_out :: Float64          # Outlet O₂ flux at the cathode GC exit

    function GCOxygenFlows{NB_GC}(cgc_in, cgc_cgc, cgc_out) where {NB_GC}
        length(cgc_cgc) == NB_GC - 1 || throw(ArgumentError("cgc_cgc length must be NB_GC - 1 = $(NB_GC - 1)."))
        return new{NB_GC}(Float64(cgc_in), collect(Float64, cgc_cgc), Float64(cgc_out))
    end
end

"""Along-channel nitrogen molar fluxes at both GC sides. Units: mol·m⁻²·s⁻¹

Fields agc_agc and cgc_cgc are inter-node flux vectors of length NB_GC - 1.
"""
struct GCNitrogenFlows{NB_GC}
    agc_in  :: Float64          # Inlet N₂ flux at the anode GC entry
    agc_agc :: Vector{Float64}  # Along-channel inter-node N₂ fluxes at anode (length NB_GC - 1)
    agc_out :: Float64          # Outlet N₂ flux at the anode GC exit
    cgc_in  :: Float64          # Inlet N₂ flux at the cathode GC entry
    cgc_cgc :: Vector{Float64}  # Along-channel inter-node N₂ fluxes at cathode (length NB_GC - 1)
    cgc_out :: Float64          # Outlet N₂ flux at the cathode GC exit

    function GCNitrogenFlows{NB_GC}(agc_in, agc_agc, agc_out,
                                     cgc_in, cgc_cgc, cgc_out) where {NB_GC}
        length(agc_agc) == NB_GC - 1 || throw(ArgumentError("agc_agc length must be NB_GC - 1 = $(NB_GC - 1)."))
        length(cgc_cgc) == NB_GC - 1 || throw(ArgumentError("cgc_cgc length must be NB_GC - 1 = $(NB_GC - 1)."))
        return new{NB_GC}(Float64(agc_in), collect(Float64, agc_agc), Float64(agc_out),
                          Float64(cgc_in), collect(Float64, cgc_cgc), Float64(cgc_out))
    end
end

"""Total (global) mass flow rates through the GC inlets and outlets. Units: mol·s⁻¹"""
struct GCMassFlows
    Wa_in  :: Float64   # Total molar flow entering the anode GC
    Wa_out :: Float64   # Total molar flow leaving the anode GC
    Wc_in  :: Float64   # Total molar flow entering the cathode GC
    Wc_out :: Float64   # Total molar flow leaving the cathode GC
end

"""Desired inlet flow setpoint at stack level. Units: mol·s⁻¹"""
struct DesiredInletFlows
    H2       :: Float64
    dry_air  :: Float64
    H2O_inj_a::Float64
    H2O_inj_c::Float64
end

"""Complete along-channel flow container for one simulation step.

Groups all GC-level typed flow structs returned by `calculate_flows_1D_GC_manifold`.
"""
struct GCManifoldFlows1D{NB_GC}
    Jv   :: GCVaporFlows{NB_GC}
    Jl   :: GCLiquidFlows{NB_GC}
    J_H2 :: GCHydrogenFlows{NB_GC}
    J_O2 :: GCOxygenFlows{NB_GC}
    J_N2 :: GCNitrogenFlows{NB_GC}
    W    :: GCMassFlows
end

"""Complete 1D MEA flow output for one gas-channel column.
"""
struct MEAFlows1D{NB_GDL, NB_MPL}
    Jv       :: MEAVaporFluxes{NB_GDL, NB_MPL}
    Jl       :: MEALiquidFluxes{NB_GDL, NB_MPL}
    J_lambda :: MEADissolvedWaterFlux
    J_H2     :: MEAHydrogenFluxes{NB_GDL, NB_MPL}
    J_O2     :: MEAOxygenFluxes{NB_GDL, NB_MPL}
    S_abs    :: MEASorptionSources
    Sp       :: MEAWaterProductionSources
    S_H2     :: MEAGasReactionSources
    S_O2     :: MEAGasReactionSources
    Sv       :: MEAVaporSources{NB_GDL, NB_MPL}
    Sl       :: MEALiquidSources{NB_GDL, NB_MPL}
end

"""Complete 1D MEA heat transfer output for one gas-channel column.

This is the future return type of `calculate_heat_transfers`, replacing
the former `Dict{String, Dict}`.
"""
struct MEAHeatFlows1D{NB_GDL, NB_MPL}
    Jt     :: MEAThermalFluxes{NB_GDL, NB_MPL}
    Q_r    :: MEAReactionHeat
    Q_sorp :: MEASorptionHeat
    Q_liq  :: MEALiquidHeat{NB_GDL, NB_MPL}
    Q_p    :: MEAProtonHeat
    Q_e    :: MEAElectricHeat{NB_GDL, NB_MPL}
end
