# -*- coding: utf-8 -*-

"""
    EpsilonBoundRestriction

Conservative pre-restriction of the catalyst-layer parameter bounds so that
`epsilon_carb`, `epsilon_mc` and `epsilon_cl` remain inside `[0, 1]` for every
sample drawn by the validity-analysis pipeline.

Only the electrode-specific parameters that enter the volume-fraction formulas
are touched:

    Hccl, IC_ccl, wt_Pt_ccl, L_Pt_ccl   (cathode)
    Hacl, IC_acl, wt_Pt_acl, L_Pt_acl   (anode, if present)

The method scales the relevant half-intervals toward their nominal stack values
by a common factor `s` chosen by bisection so that the worst-case total solid +
ionomer fraction is exactly `1`.  This yields a conservative rectangular box
inscribed in the feasible region; it is safe but not volume-optimal.
"""
module EpsilonBoundRestriction

using Printf
using AlphaPEM.Config: PhysicalParams
using AlphaPEM.Fuelcell: create_fuelcell, undetermined_parameter_bounds
using AlphaPEM.Utils: rho_H2O_l, M_H2O, rho_carb, rho_Pt

export restrict_epsilon_bounds, verify_restricted_bounds

"""Return `(element, property)` for the geometrical parameters that influence
the catalyst-layer volume fractions.  `element` is `:acl` or `:ccl`;
`property` is `:Hcl`, `:IC`, `:wt` or `:L_Pt`."""
function parameter_role(name::Symbol)
    name == :Hccl       && return (:ccl, :Hcl)
    name == :Hacl       && return (:acl, :Hcl)
    name == :IC_ccl     && return (:ccl, :IC)
    name == :IC_acl     && return (:acl, :IC)
    name == :wt_Pt_ccl  && return (:ccl, :wt)
    name == :wt_Pt_acl  && return (:acl, :wt)
    name == :L_Pt_ccl   && return (:ccl, :L_Pt)
    name == :L_Pt_acl   && return (:acl, :L_Pt)
    return nothing
end

function element_param_name(element::Symbol, property::Symbol)::Symbol
    if element == :ccl
        property == :Hcl  && return :Hccl
        property == :IC   && return :IC_ccl
        property == :wt   && return :wt_Pt_ccl
        property == :L_Pt && return :L_Pt_ccl
    else
        property == :Hcl  && return :Hacl
        property == :IC   && return :IC_acl
        property == :wt   && return :wt_Pt_acl
        property == :L_Pt && return :L_Pt_acl
    end
    throw(ArgumentError("Unknown (element,property) = ($element,$property)"))
end

"""Worst-case value of the ionomer swelling factor K(T,lambda)."""
function max_K(pp::PhysicalParams, lambda_range::Tuple{Float64,Float64},
               T_range::Tuple{Float64,Float64})::Float64
    lambda_max = max(lambda_range...)
    T_max      = max(T_range...)
    return 1.0 + (M_H2O * pp.rho_ion) / (rho_H2O_l(T_max) * pp.M_eq) * lambda_max
end

"""Worst-case total fraction over a rectangular parameter box."""
function worst_case_total(element::Symbol, pp::PhysicalParams, Kmax::Float64,
                          bounds::Dict{Symbol,Tuple{Float64,Float64}})::Float64
    function getb(element, prop)
        name = element_param_name(element, prop)
        if haskey(bounds, name)
            return bounds[name]
        end
        val = getproperty(pp, name)
        return (val, val)
    end

    Hcl_min, _  = getb(element, :Hcl)
    _, IC_max   = getb(element, :IC)
    wt_min, _   = getb(element, :wt)
    _, L_Pt_max = getb(element, :L_Pt)

    return L_Pt_max / Hcl_min * (
        (1 - wt_min) / (wt_min * rho_carb) +
        1 / rho_Pt +
        IC_max * (1 - wt_min) / (wt_min * pp.rho_ion) * Kmax
    )
end

"""
    restrict_epsilon_bounds(pp, original_bounds;
                            lambda_range=(0.0, 22.0),
                            T_range=(273.15 + 20.0, 273.15 + 95.0),
                            verbose=false)

Return a new dictionary of conservative bounds such that, for every parameter
sample drawn inside the new box and for every `(lambda,T)` inside the given
ranges, the fractions `epsilon_carb`, `epsilon_mc` and `epsilon_cl` stay in
`[0,1]` for both catalyst layers.
"""
function restrict_epsilon_bounds(pp::PhysicalParams,
                                 original_bounds::Dict{Symbol,Tuple{Float64,Float64}};
                                 lambda_range::Tuple{Float64,Float64}=(0.0, 22.0),
                                 T_range::Tuple{Float64,Float64}=(273.15 + 20.0,
                                                                   273.15 + 95.0),
                                 verbose::Bool=false)

    Kmax = max_K(pp, lambda_range, T_range)
    new_bounds = copy(original_bounds)

    electrodes = Set{Symbol}()
    for name in keys(original_bounds)
        role = parameter_role(name)
        isnothing(role) && continue
        push!(electrodes, role[1])
    end
    push!(electrodes, :acl)
    push!(electrodes, :ccl)

    for element in electrodes
        S_wc = worst_case_total(element, pp, Kmax, new_bounds)

        if verbose
            elname = element == :ccl ? "CCL" : "ACL"
            @printf("  %s worst-case total fraction = %.5f (must be <= 1.0)\n",
                    elname, S_wc)
        end

        S_wc <= 1.0 && continue

        relevant = [p for p in (:Hcl, :IC, :wt, :L_Pt)
                      if haskey(new_bounds, element_param_name(element, p))]
        if isempty(relevant)
            error("Known parameters for $element already violate the epsilon " *
                  "constraint (S=$S_wc) and there is no uncertain parameter to " *
                  "restrict. Check the stack data.")
        end

        nominal_val(prop) = getproperty(pp, element_param_name(element, prop))
        function orig_bounds(prop)
            n = element_param_name(element, prop)
            haskey(new_bounds, n) ? new_bounds[n] : (nominal_val(prop), nominal_val(prop))
        end

        Hcl_nom = nominal_val(:Hcl)
        IC_nom  = nominal_val(:IC)
        wt_nom  = nominal_val(:wt)
        L_nom   = nominal_val(:L_Pt)

        Hcl_lo, Hcl_hi = orig_bounds(:Hcl)
        IC_lo,  IC_hi  = orig_bounds(:IC)
        wt_lo,  wt_hi  = orig_bounds(:wt)
        L_lo,   L_hi   = orig_bounds(:L_Pt)

        function S_of_s(s::Float64)::Float64
            Lmax  = L_hi   > L_nom   ? L_nom   + s * (L_hi   - L_nom)   : L_hi
            ICmax = IC_hi  > IC_nom  ? IC_nom  + s * (IC_hi  - IC_nom)  : IC_hi
            Hmin  = Hcl_lo < Hcl_nom ? Hcl_nom - s * (Hcl_nom - Hcl_lo) : Hcl_lo
            wtmin = wt_lo  < wt_nom  ? wt_nom  - s * (wt_nom  - wt_lo)  : wt_lo

            tmp = copy(new_bounds)
            tmp[element_param_name(element, :L_Pt)] = (L_lo,  Lmax)
            tmp[element_param_name(element, :IC)]   = (IC_lo, ICmax)
            tmp[element_param_name(element, :Hcl)]  = (Hmin,  Hcl_hi)
            tmp[element_param_name(element, :wt)]   = (wtmin, wt_hi)
            return worst_case_total(element, pp, Kmax, tmp)
        end

        S0 = S_of_s(0.0)
        if S0 > 1.0
            error("Nominal parameters already violate the epsilon constraint for $element (S=$S0).")
        end

        s_lo, s_hi = 0.0, 1.0
        for _ in 1:60
            s_mid = (s_lo + s_hi) / 2
            if S_of_s(s_mid) <= 1.0
                s_lo = s_mid
            else
                s_hi = s_mid
            end
        end
        s = s_lo

        Lmax  = L_hi   > L_nom   ? L_nom   + s * (L_hi   - L_nom)   : L_hi
        ICmax = IC_hi  > IC_nom  ? IC_nom  + s * (IC_hi  - IC_nom)  : IC_hi
        Hmin  = Hcl_lo < Hcl_nom ? Hcl_nom - s * (Hcl_nom - Hcl_lo) : Hcl_lo
        wtmin = wt_lo  < wt_nom  ? wt_nom  - s * (wt_nom  - wt_lo)  : wt_lo

        name_L  = element_param_name(element, :L_Pt)
        name_IC = element_param_name(element, :IC)
        name_H  = element_param_name(element, :Hcl)
        name_wt = element_param_name(element, :wt)
        haskey(new_bounds, name_L)  && (new_bounds[name_L]  = (L_lo,  Lmax))
        haskey(new_bounds, name_IC) && (new_bounds[name_IC] = (IC_lo, ICmax))
        haskey(new_bounds, name_H)  && (new_bounds[name_H]  = (Hmin,  Hcl_hi))
        haskey(new_bounds, name_wt) && (new_bounds[name_wt] = (wtmin, wt_hi))

        S_new = worst_case_total(element, pp, Kmax, new_bounds)
        if S_new > 1.0 + 1e-12
            error("Internal error: restriction failed for $element (S=$S_new).")
        end
    end

    return new_bounds
end

"""
    verify_restricted_bounds(pp, bounds, lambda_range, T_range; n_samples=2000)

Draw `n_samples` random combinations of parameters, lambda and T inside the
boxes and ranges, recompute the volume fractions with the original AlphaPEM
formulas, and check that they all lie in `[0,1]`.
"""
function verify_restricted_bounds(pp::PhysicalParams,
                                  bounds::Dict{Symbol,Tuple{Float64,Float64}},
                                  lambda_range::Tuple{Float64,Float64},
                                  T_range::Tuple{Float64,Float64};
                                  n_samples::Int=2000)::Bool

    function sample_box()
        d = Dict{Symbol,Float64}()
        for (name, (lo, hi)) in bounds
            d[name] = lo + rand() * (hi - lo)
        end
        return d
    end

    for _ in 1:n_samples
        s = sample_box()
        lambda = lambda_range[1] + rand() * (lambda_range[2] - lambda_range[1])
        T      = T_range[1]      + rand() * (T_range[2]      - T_range[1])

        for element in (:acl, :ccl)
            function val(prop)
                n = element_param_name(element, prop)
                haskey(s, n) && return s[n]
                return getproperty(pp, n)
            end
            Hcl, IC, wt, L_Pt = val(:Hcl), val(:IC), val(:wt), val(:L_Pt)

            eps_carb = L_Pt * (1 - wt) / wt / (rho_carb * Hcl)
            eps_Pt   = L_Pt / (rho_Pt * Hcl)
            K = 1.0 + (M_H2O * pp.rho_ion) / (rho_H2O_l(T) * pp.M_eq) * lambda
            eps_mc   = IC * eps_carb * rho_carb / pp.rho_ion * K
            eps_cl   = 1.0 - eps_carb - eps_Pt - eps_mc

            if !(0.0 <= eps_carb <= 1.0) || !(0.0 <= eps_mc <= 1.0) || !(0.0 <= eps_cl <= 1.0)
                @warn("Violation for $element: " *
                      "eps_carb=$eps_carb, eps_mc=$eps_mc, eps_cl=$eps_cl\n" *
                      "sample: Hcl=$Hcl, IC=$IC, wt=$wt, L_Pt=$L_Pt, T=$T, lambda=$lambda")
                return false
            end
        end
    end
    return true
end

end  # module EpsilonBoundRestriction
