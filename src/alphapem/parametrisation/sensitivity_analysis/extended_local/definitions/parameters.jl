# -*- coding: utf-8 -*-

"""
    PARAMETER_LATEX

Display symbol (LaTeX) of each parameter of `PARAMETER_BOUNDS`. Parameters without an entry fall
back to `string(name)`. Provisional table, to be merged into the central parameter table.
"""
const PARAMETER_LATEX = Dict{Symbol, String}(
    :Hacl => "H_{acl}", :Hccl => "H_{ccl}", :Hmem => "H_{mem}", :Hgdl => "H_{gdl}",
    :epsilon_gdl => "\\varepsilon_{gdl}", :Hmpl => "H_{mpl}", :epsilon_mpl => "\\varepsilon_{mpl}",
    :Hagc => "H_{agc}", :Hcgc => "H_{cgc}", :Wagc => "W_{agc}", :Wcgc => "W_{cgc}", :Lgc => "L_{gc}",
    :rho_mem => "\\rho_{mem}", :M_eq => "M_{eq}", :tau_mpl => "\\tau_{mpl}", :tau_cl => "\\tau_{cl}",
    :r_s_gdl => "r_{s,gdl}", :r_s_mpl => "r_{s,mpl}", :r_s_cl => "r_{s,cl}",
    :Dp_gdl => "D_{p,gdl}", :Dp_mpl => "D_{p,mpl}", :Dp_cl => "D_{p,cl}",
    :theta_c_gdl => "\\theta_{c,gdl}", :theta_c_mpl => "\\theta_{c,mpl}", :theta_c_cl => "\\theta_{c,cl}",
    :theta_l_rem => "\\theta_{l,rem}", :K_v_liq_gas => "K_{v,l/g}", :gamma_sorp_l => "\\gamma_{sorp,l}",
    :gamma_cond => "\\gamma_{cond}", :gamma_evap => "\\gamma_{evap}",
    :epsilon_p => "\\varepsilon_{p}", :epsilon_c => "\\varepsilon_{c}", :alpha_p => "\\alpha_{p}",
    :Eact_H2_cros_v => "E_{act,H_2,cros,v}", :Eact_H2_cros_l => "E_{act,H_2,cros,l}",
    :Eact_O2_cros_v => "E_{act,O_2,cros,v}", :Eact_O2_cros_l => "E_{act,O_2,cros,l}",
    :e => "e", :K_O2_ad_Pt => "K_{O_2,ad,Pt}", :K_O2_dis_ion => "K_{O_2,dis,ion}",
    :K_O2_dis_l => "K_{O_2,dis,l}", :IC_acl => "IC_{acl}", :IC_ccl => "IC_{ccl}",
    :ECSA_0 => "ECSA_0", :wt_Pt_acl => "wt_{Pt,acl}", :wt_Pt_ccl => "wt_{Pt,ccl}",
    :L_Pt_acl => "L_{Pt,acl}", :L_Pt_ccl => "L_{Pt,ccl}", :r_carb => "r_{carb}", :rho_ion => "\\rho_{ion}",
    :Re => "R_e", :i0_c_ref => "i_{0,c,ref}", :kappa_co => "\\kappa_{co}", :kappa_c => "\\kappa_{c}",
    :Eact_O2_red => "E_{act,O_2,red}", :C_scl => "C_{scl}", :alpha_c => "\\alpha_{c}",
    :k_th_gdl => "k_{th,gdl}", :k_th_mpl => "k_{th,mpl}", :k_th_cl => "k_{th,cl}", :k_th_mem => "k_{th,mem}",
    :Cp_gdl => "C_{p,gdl}", :Cp_mpl => "C_{p,mpl}", :Cp_cl => "C_{p,cl}", :Cp_mem => "C_{p,mem}",
    :rho_gdl => "\\rho_{gdl}", :rho_mpl => "\\rho_{mpl}", :rho_cl => "\\rho_{cl}",
    :sigma_e_gdl => "\\sigma_{e,gdl}", :sigma_e_mpl => "\\sigma_{e,mpl}", :sigma_e_cl => "\\sigma_{e,cl}",
)


"""
    all_parameter_bounds() -> Vector{ParameterBound}

Build parameter bounds and metadata for every entry of `PARAMETER_BOUNDS`.

Parameters
----------
None

Returns
-------
Vector{ParameterBound}
    Bounds sorted by parameter name, with scale and display metadata.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function all_parameter_bounds()::Vector{ParameterBound}
    # Prepare the sorted result container
    bounds = ParameterBound[]

    # Validate each bound and combine it with scale and display metadata
    for name in sort!(collect(keys(PARAMETER_BOUNDS)))
        lo, hi, ptype = PARAMETER_BOUNDS[name]
        haskey(PARAMETER_METADATA, name) ||
            throw(ArgumentError("No metadata found for parameter $name"))
        lo < hi || throw(ArgumentError("Invalid bounds for $name: min=$lo >= max=$hi"))
        unit, description = PARAMETER_METADATA[name]
        scale = Calibration._default_parameter_scale(Float64(lo), Float64(hi), ptype)
        push!(bounds, ParameterBound(name, Float64(lo), Float64(hi), ptype, scale, unit, description,
                                     get(PARAMETER_LATEX, name, string(name))))
    end

    # Return the complete parameter list
    return bounds
end


"""
    select_parameters(cfg::ExtendedLocalSensitivityConfig) -> (params, fixed)

Select the parameters to analyse and identify integer parameters held fixed.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Analysis configuration, including optional bounds, subset, and integer policy.

Returns
-------
Tuple{Vector{ParameterBound}, Vector{Symbol}}
    The selected parameter bounds and names of integer parameters kept at nominal values.

Raises
------
ArgumentError
    If the requested subset contains unknown parameters or selects no analysable parameters.
"""
function select_parameters(cfg::ExtendedLocalSensitivityConfig)
    # Start with custom bounds or the complete central parameter list
    params = cfg.parameter_bounds === nothing ? all_parameter_bounds() : copy(cfg.parameter_bounds.bounds)
    params = [p.latex == "" ? ParameterBound(p.name, p.min, p.max, p.type, p.scale, p.unit, p.description,
                                             get(PARAMETER_LATEX, p.name, string(p.name))) : p
              for p in params]

    # Restrict the list to the requested subset
    if cfg.parameter_subset !== nothing
        known = Set(p.name for p in params)
        unknown = [n for n in cfg.parameter_subset if !(n in known)]
        isempty(unknown) || throw(ArgumentError("Unknown parameters in parameter_subset: $unknown"))
        keep = Set(cfg.parameter_subset)
        params = filter(p -> p.name in keep, params)
    end

    # Apply the configured policy for integer parameters
    fixed = Symbol[]
    if cfg.integer_policy == :fix
        fixed = [p.name for p in params if p.type == :int]
        isempty(fixed) || @warn "Integer parameters kept at their nominal value and excluded: $fixed"
        params = filter(p -> p.type != :int, params)
    end

    # Ensure the analysis has at least one free parameter
    isempty(params) && throw(ArgumentError("No parameter left to analyse"))
    return params, fixed
end


"""
    to_unit(p::ParameterBound, θ)

Map physical parameter values to the unit interval using the parameter's scale.

Parameters
----------
p : ParameterBound
    Parameter bounds and linear/logarithmic scale.
θ : Real or array-like
    Physical value or values to map.

Returns
-------
Float64 or array
    Value or values mapped to `[0, 1]`.

Raises
------
ArgumentError
    If logarithmic scaling is requested for a parameter with a non-positive lower bound.
"""
function to_unit(p::ParameterBound, θ)
    # Apply the logarithmic mapping when the parameter uses a log scale
    if p.scale == :log
        p.min > 0.0 || throw(ArgumentError("Log scale requires min > 0 for $(p.name) (got $(p.min))"))
        return @. (log(θ) - log(p.min)) / (log(p.max) - log(p.min))
    end

    # Otherwise use a linear mapping
    return @. (θ - p.min) / (p.max - p.min)
end


"""
    from_unit(p::ParameterBound, u; round_int::Bool = false)

Map unit-interval coordinates to physical parameter values.

Parameters
----------
p : ParameterBound
    Parameter bounds and linear/logarithmic scale.
u : Real or array-like
    Unit-interval coordinate or coordinates.
round_int : Bool, optional
    Round integer parameters to the nearest integer. Default is `false`.

Returns
-------
Float64 or array
    Physical value or values corresponding to `u`.

Raises
------
ArgumentError
    If logarithmic scaling is requested for a parameter with a non-positive lower bound.
"""
function from_unit(p::ParameterBound, u; round_int::Bool = false)
    # Map unit coordinates back through the parameter's selected scale
    if p.scale == :log
        p.min > 0.0 || throw(ArgumentError("Log scale requires min > 0 for $(p.name) (got $(p.min))"))
        θ = @. exp(log(p.min) + u * (log(p.max) - log(p.min)))
    else
        θ = @. p.min + u * (p.max - p.min)
    end

    # Round only when requested for an integer-valued parameter
    return (round_int && p.type == :int) ? round.(θ) : θ
end


"""
    _physical_params_with_overrides(base, params, values) -> PhysicalParams

Return `base` with the analysed parameters replaced by `values`, all other fields unchanged.

Parameters
----------
base : PhysicalParams
    Nominal physical parameters.
params : AbstractVector
    Analysed parameters, each with the fields `name` and `type`.
values : AbstractVector{<:Real}
    Physical values in the order of `params`; integer parameters are rounded.

Returns
-------
PhysicalParams
    New parameter set.
"""
function _physical_params_with_overrides(base::PhysicalParams, params::AbstractVector,
                                         values::AbstractVector{<:Real})::PhysicalParams
    base_fields = (; (f => getfield(base, f) for f in fieldnames(PhysicalParams))...)
    overrides = (; (p.name => (p.type == :int ? Int(round(v)) : Float64(v))
                    for (p, v) in zip(params, values))...)
    return PhysicalParams(; merge(base_fields, overrides)...)
end
