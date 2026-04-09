# This file represents the differential equations for the gas channels and manifolds.
#
# Each public function computes one physics contribution (gas species, liquid water, temperature)
# and returns per-node Float64 vectors. The caller assembles AnodeGCDerivative /
# CathodeGCDerivative and then the full MEACellDerivative1D.

# ____________________________________________________Main functions____________________________________________________

"""Assemble complete MEA derivatives by injecting GC contributions into existing MEA derivatives."""
function assemble_gc_derivative_1D(mea::AbstractVector{<:MEACellDerivative1D{NB_GDL, NB_MPL}},
                                   gas::GCGasDerivative{NB_GC},
                                   liq::GCLiquidWaterDerivative{NB_GC},
                                   temp::GCTemperatureDerivative{NB_GC}) where {NB_GDL, NB_MPL, NB_GC}
    return [begin
                agc = AnodeGCDerivative(gas.agc_C_v[i], liq.agc_s[i], gas.agc_C_H2[i], gas.agc_C_N2[i], temp.agc_T[i])
                cgc = CathodeGCDerivative(gas.cgc_C_v[i], liq.cgc_s[i], gas.cgc_C_O2[i], gas.cgc_C_N2[i], temp.cgc_T[i])
                node = mea[i]
                MEACellDerivative1D{NB_GDL, NB_MPL}(agc, node.agdl, node.ampl, node.acl, node.mem,
                                                    node.ccl, node.cmpl, node.cgdl, cgc)
            end for i in 1:NB_GC]
end

"""Calculate dynamic gas-species (vapour, H₂, O₂, N₂) evolution in the gas channels.

Returns per-node derivative vectors for each species at the anode and cathode GC.
GC-to-MEA interface fluxes are read directly from `flows_mea` (MEAFlows1D per node).

Parameters
----------
sv             : typed state vector (one MEAState1D per GC node)
pp             : physical-parameters container (`pp.Hagc`, `pp.Hcgc`, `pp.Lgc`)
cfg            : simulation configuration (`cfg.type_auxiliary`)
flows_gc       : GCManifoldFlows1D — along-channel GC fluxes
flows_mea      : Vector of MEAFlows1D — one entry per GC node (GC-MEA interface fluxes)

Returns
-------
GCGasDerivative{nb_gc}
    Typed gas-species derivative contribution for the gas channels.
"""
function calculate_dyn_gas_evolution_inside_gas_channel(
        sv             :: AbstractVector{<:MEAState1D},
        pp             :: PhysicalParams,
        cfg            :: SimulationConfig,
        flows_gc       :: GCManifoldFlows1D{NB_GC},
        flows_mea      :: AbstractVector) where {NB_GC}

    Jv  = flows_gc.Jv
    JH2 = flows_gc.J_H2
    JO2 = flows_gc.J_O2
    JN2 = flows_gc.J_N2
    Hagc = pp.Hagc
    Hcgc = pp.Hcgc
    Lgc = pp.Lgc
    type_auxiliary = cfg.type_auxiliary
    L_node = Lgc / NB_GC

    # Anode GC: water vapour
    agc_C_v = ntuple(NB_GC) do i
        fac_a = 1.0 / (1.0 - sv[i].agc.s)
        if NB_GC == 1
            J_in = Jv.agc_in
            J_out = Jv.agc_out
            fac_a * (J_in - J_out) / Lgc - flows_mea[i].Jv.agc_agdl / Hagc
        else
            J_in = i == 1 ? Jv.agc_in : Jv.agc_agc[i - 1]
            J_out = i == NB_GC ? Jv.agc_out : Jv.agc_agc[i]
            fac_a * (J_in - J_out) / L_node - flows_mea[i].Jv.agc_agdl / Hagc
        end
    end

    # Anode GC: hydrogen
    agc_C_H2 = ntuple(NB_GC) do i
        fac_a = 1.0 / (1.0 - sv[i].agc.s)
        if NB_GC == 1
            J_in = JH2.agc_in
            J_out = JH2.agc_out
            fac_a * (J_in - J_out) / Lgc - flows_mea[i].J_H2.agc_agdl / Hagc
        else
            J_in = i == 1 ? JH2.agc_in : JH2.agc_agc[i - 1]
            J_out = i == NB_GC ? JH2.agc_out : JH2.agc_agc[i]
            fac_a * (J_in - J_out) / L_node - flows_mea[i].J_H2.agc_agdl / Hagc
        end
    end

    # Anode GC: nitrogen
    agc_C_N2 = ntuple(NB_GC) do i
        fac_a = 1.0 / (1.0 - sv[i].agc.s)
        if type_auxiliary != :forced_convective_cathode_with_flow_through_anode
            0.0
        elseif NB_GC == 1
            J_in = JN2.agc_in
            J_out = JN2.agc_out
            fac_a * (J_in - J_out) / Lgc
        else
            J_in = i == 1 ? JN2.agc_in : JN2.agc_agc[i - 1]
            J_out = i == NB_GC ? JN2.agc_out : JN2.agc_agc[i]
            fac_a * (J_in - J_out) / L_node
        end
    end

    # Cathode GC: water vapour
    cgc_C_v = ntuple(NB_GC) do i
        fac_c = 1.0 / (1.0 - sv[i].cgc.s)
        if NB_GC == 1
            J_in = Jv.cgc_in
            J_out = Jv.cgc_out
            fac_c * (J_in - J_out) / Lgc + flows_mea[i].Jv.cgdl_cgc / Hcgc
        else
            J_in = i == 1 ? Jv.cgc_in : Jv.cgc_cgc[i - 1]
            J_out = i == NB_GC ? Jv.cgc_out : Jv.cgc_cgc[i]
            fac_c * (J_in - J_out) / L_node + flows_mea[i].Jv.cgdl_cgc / Hcgc
        end
    end

    # Cathode GC: oxygen
    cgc_C_O2 = ntuple(NB_GC) do i
        fac_c = 1.0 / (1.0 - sv[i].cgc.s)
        if NB_GC == 1
            J_in = JO2.cgc_in
            J_out = JO2.cgc_out
            fac_c * (J_in - J_out) / Lgc + flows_mea[i].J_O2.cgdl_cgc / Hcgc
        else
            J_in = i == 1 ? JO2.cgc_in : JO2.cgc_cgc[i - 1]
            J_out = i == NB_GC ? JO2.cgc_out : JO2.cgc_cgc[i]
            fac_c * (J_in - J_out) / L_node + flows_mea[i].J_O2.cgdl_cgc / Hcgc
        end
    end

    # Cathode GC: nitrogen
    cgc_C_N2 = ntuple(NB_GC) do i
        fac_c = 1.0 / (1.0 - sv[i].cgc.s)
        if NB_GC == 1
            J_in = JN2.cgc_in
            J_out = JN2.cgc_out
            fac_c * (J_in - J_out) / Lgc
        else
            J_in = i == 1 ? JN2.cgc_in : JN2.cgc_cgc[i - 1]
            J_out = i == NB_GC ? JN2.cgc_out : JN2.cgc_cgc[i]
            fac_c * (J_in - J_out) / L_node
        end
    end

    return GCGasDerivative{NB_GC}(agc_C_v, agc_C_H2, agc_C_N2, cgc_C_v, cgc_C_O2, cgc_C_N2)
end


"""Calculate dynamic liquid-water evolution in the gas channels.

Returns per-node derivative vectors for the liquid saturation at anode and cathode GC.

Parameters
----------
T_des       : desired fuel-cell temperature (K)
pp          : physical-parameters container (`pp.Hagc`, `pp.Hcgc`, `pp.Lgc`)
flows_gc    : GCManifoldFlows1D — along-channel liquid GC fluxes
flows_mea   : Vector of MEAFlows1D — one entry per GC node (GC-MEA interface fluxes)

Returns
-------
GCLiquidWaterDerivative{nb_gc}
    Typed liquid-water derivative contribution for the gas channels.
"""
function calculate_dyn_liq_evolution_inside_gas_channel(
        T_des     :: Float64,
        pp        :: PhysicalParams,
        flows_gc  :: GCManifoldFlows1D{NB_GC},
        flows_mea :: AbstractVector) where {NB_GC}

    inv_rho = 1.0 / rho_H2O_l(T_des)
    Jl = flows_gc.Jl
    Hagc = pp.Hagc
    Hcgc = pp.Hcgc
    Lgc = pp.Lgc
    L_node = Lgc / NB_GC

    # Anode GC: liquid water
    agc_s = ntuple(NB_GC) do i
        if NB_GC == 1
            J_in = 0.0
            J_out = Jl.agc_out
            inv_rho * ((J_in - J_out) / Lgc - flows_mea[i].Jl.agc_agdl / Hagc)
        else
            J_in = i == 1 ? 0.0 : Jl.agc_agc[i - 1]
            J_out = i == NB_GC ? Jl.agc_out : Jl.agc_agc[i]
            inv_rho * ((J_in - J_out) / L_node - flows_mea[i].Jl.agc_agdl / Hagc)
        end
    end

    # Cathode GC: liquid water
    cgc_s = ntuple(NB_GC) do i
        if NB_GC == 1
            J_in = 0.0
            J_out = Jl.cgc_out
            inv_rho * ((J_in - J_out) / Lgc + flows_mea[i].Jl.cgdl_cgc / Hcgc)
        else
            J_in = i == 1 ? 0.0 : Jl.cgc_cgc[i - 1]
            J_out = i == NB_GC ? Jl.cgc_out : Jl.cgc_cgc[i]
            inv_rho * ((J_in - J_out) / L_node + flows_mea[i].Jl.cgdl_cgc / Hcgc)
        end
    end

    return GCLiquidWaterDerivative{NB_GC}(agc_s, cgc_s)
end


"""Calculate dynamic temperature evolution in the gas channels.

Temperature in both GCs is imposed by a Dirichlet boundary condition (initialized to T_fc,
kept constant), so all derivatives are identically zero.

Parameters
----------
nb_gc : number of gas channels (integer)

Returns
-------
GCTemperatureDerivative{nb_gc}
    Typed temperature derivative contribution for the gas channels.
"""
function calculate_dyn_temperature_evolution_inside_gas_channel(nb_gc)
    return GCTemperatureDerivative{nb_gc}(ntuple(_ -> 0.0, nb_gc), ntuple(_ -> 0.0, nb_gc))
end


"""Calculate pressure and humidity dynamics in the manifolds.

Parameters
----------
dif_eq_manifold_1D : typed manifold derivative container (mutable or NaN placeholder)
T_des              : fuel cell temperature (K)
type_auxiliary     : auxiliary system type (Symbol)
W                  : GCMassFlows — total molar inlet/outlet flows (mol·s⁻¹)
Wv                 : vapour molar flows for the manifolds (not yet fully implemented)
"""
function calculate_dyn_manifold_pressure_and_humidity_evolution(
    _dif_eq_manifold_1D,
    _T_des          :: Float64,
    _type_auxiliary :: Symbol,
    _W              :: GCMassFlows,
    _Wv
)

    # # Pressure evolution inside the manifolds
    # if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
    #    type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #     # At the anode side
    #     if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
    #         dif_eq["dPasm / dt"] = (W.Wa_in + W_asm_in_re_to_asm - nb_cell * Wasm_to_asm_out) / Vasm * R * T_des
    #         dif_eq["dPaem / dt"] = (nb_cell * Waem_in_to_aem - Waem_to_aem_out - Waem_to_aem_out_re) / Vaem * R * T_des
    #     else  # type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #         dif_eq["dPasm / dt"] = (W.Wa_in - nb_cell * Wasm_to_asm_out) / Vasm * R * T_des
    #         dif_eq["dPaem / dt"] = (nb_cell * Waem_in_to_aem - Waem_to_aem_out) / Vaem * R * T_des
    #     end
    #     # At the cathode side
    #     dif_eq["dPcsm / dt"] = (W.Wc_in - nb_cell * Wcsm_to_csm_out) / Vcsm * R * T_des
    #     dif_eq["dPcem / dt"] = (nb_cell * Wcem_in_to_cem - Wcem_to_cem_out) / Vcem * R * T_des
    # end
    #
    # # Humidity evolution inside the manifolds
    # if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
    #    type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #     # At the anode side
    #     if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
    #         dif_eq["dPhi_asm / dt"] = (Wv.Wva_in + Wv_asm_in_re_to_asm - nb_cell * Wv_asm_to_asm_out) / Vasm * R * T_des / Psat(T_des)
    #         dif_eq["dPhi_aem / dt"] = (nb_cell * Wv_aem_in_to_aem - Wv_aem_to_aem_out_re - Wv_aem_to_aem_out) / Vaem * R * T_des / Psat(T_des)
    #     else  # type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #         dif_eq["dPhi_asm / dt"] = (Wv.Wva_in - nb_cell * Wv_asm_to_asm_out) / Vasm * R * T_des / Psat(T_des)
    #         dif_eq["dPhi_aem / dt"] = (nb_cell * Wv_aem_in_to_aem - Wv_aem_to_aem_out) / Vaem * R * T_des / Psat(T_des)
    #     end
    #     # At the cathode side
    #     dif_eq["dPhi_csm / dt"] = (Wv.Wvc_in - nb_cell * Wv_csm_to_csm_out) / Vcsm * R * T_des / Psat(T_des)
    #     dif_eq["dPhi_cem / dt"] = ...
    # end
    return nothing
end

