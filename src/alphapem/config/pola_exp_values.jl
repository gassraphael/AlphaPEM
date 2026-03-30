# -*- coding: utf-8 -*-

"""This file is designated for experimental polarization-curve data used by AlphaPEM."""

# ___________________________________________________Experimental data__________________________________________________

"""
    pola_exp_values_calibration(type_fuel_cell, voltage_zone)

This function returns the experimental values of polarisation curves made on different fuel cells at different
operating conditions. The experimental values are used for calibrating the model and so are composed of a reduced
number of points compare to the pola_exp_values function. These points are specifically chosen to be as few as
possible while still providing a good representation of the polarisation curve.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell used in the model. This parameter includes the fuel cell used
  in the model and the corresponding operating conditions.
- `voltage_zone::String`: Zone of the polarization curve which is considered. It can be `"full"` or
  `"before_voltage_drop"`.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: Experimental values for calibration (current density and voltage).
"""
function pola_exp_values_calibration(type_fuel_cell::String,
                                     voltage_zone::String)::Tuple{Vector{Float64}, Vector{Float64}}

    # ZSW fuel cell
    if type_fuel_cell == "ZSW-GenStack"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099, 1.700, 2.000, 2.500]
            U_exp_cali_t = [0.953, 0.864, 0.743, 0.685, 0.620, 0.579, 0.486]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099, 1.700]
            U_exp_cali_t = [0.953, 0.864, 0.743, 0.685, 0.620]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.300, 0.700, 0.900, 1.500, 2.171]
            U_exp_cali_t = [0.936, 0.835, 0.759, 0.701, 0.670, 0.541, 0.402]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.300, 0.700]
            U_exp_cali_t = [0.936, 0.835, 0.759, 0.701]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.300, 2.000, 2.415]
            U_exp_cali_t = [0.946, 0.855, 0.736, 0.655, 0.545, 0.450]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.300]
            U_exp_cali_t = [0.946, 0.855, 0.736, 0.655]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099, 1.900, 2.200, 2.500]
            U_exp_cali_t = [0.949, 0.867, 0.746, 0.687, 0.607, 0.566, 0.514]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099, 1.900]
            U_exp_cali_t = [0.949, 0.867, 0.746, 0.687, 0.607]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099, 1.900, 2.200, 2.500]
            U_exp_cali_t = [0.947, 0.872, 0.752, 0.694, 0.622, 0.588, 0.547]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099, 1.900]
            U_exp_cali_t = [0.947, 0.872, 0.752, 0.703, 0.622]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_T_62"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.500, 2.000, 2.500]
            U_exp_cali_t = [0.944, 0.855, 0.739, 0.631, 0.566, 0.471]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.500]
            U_exp_cali_t = [0.944, 0.855, 0.739, 0.631]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_T_76"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099, 1.700, 2.200, 2.500]
            U_exp_cali_t = [0.946, 0.849, 0.744, 0.679, 0.587, 0.482, 0.406]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.498, 1.099]
            U_exp_cali_t = [0.946, 0.849, 0.744, 0.679]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_T_84"
        if voltage_zone == "full"
            i_exp_cali_t = [0.001, 0.050, 0.300, 0.700, 0.901, 1.500, 2.000]
            U_exp_cali_t = [0.930, 0.847, 0.772, 0.718, 0.686, 0.575, 0.461]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.001, 0.050, 0.300, 0.700]
            U_exp_cali_t = [0.930, 0.847, 0.772, 0.718]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    # EH-31 fuel cell
    elseif type_fuel_cell == "EH-31_1.5"  # at 1.5 bar
        if voltage_zone == "full"
            i_exp_cali_t = [0.050, 0.110, 0.293, 1.039, 1.683, 1.966, 2.246]
            U_exp_cali_t = [0.900, 0.850, 0.794, 0.681, 0.599, 0.556, 0.500]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.050, 0.110, 0.293, 1.039, 1.683]
            U_exp_cali_t = [0.900, 0.850, 0.794, 0.681, 0.599]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "EH-31_2.0"  # at 2.0 bar
        if voltage_zone == "full"
            i_exp_cali_t = [0.050, 0.106, 0.242, 0.681, 1.242, 1.501, 1.979, 2.459]
            U_exp_cali_t = [0.900, 0.860, 0.830, 0.759, 0.698, 0.668, 0.598, 0.502]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.050, 0.106, 0.242, 0.681, 1.242]
            U_exp_cali_t = [0.900, 0.860, 0.830, 0.759, 0.698]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "EH-31_2.25"  # at 2.25 bar
        if voltage_zone == "full"
            i_exp_cali_t = [0.056, 0.183, 0.364, 1.011, 1.675, 1.918, 2.356, 2.794]
            U_exp_cali_t = [0.894, 0.854, 0.819, 0.746, 0.678, 0.647, 0.580, 0.497]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.056, 0.183, 0.364, 1.011, 1.675]
            U_exp_cali_t = [0.894, 0.854, 0.819, 0.746, 0.678]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "EH-31_2.5"  # at 2.5 bar
        if voltage_zone == "full"
            i_exp_cali_t = [0.057, 0.127, 0.296, 0.527, 1.030, 1.514, 1.972, 2.358, 2.691, 2.988]
            U_exp_cali_t = [0.900, 0.866, 0.835, 0.800, 0.758, 0.712, 0.663, 0.621, 0.575, 0.524]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_cali_t = [0.057, 0.127, 0.296, 0.527, 1.030, 1.514]
            U_exp_cali_t = [0.900, 0.866, 0.835, 0.800, 0.758, 0.712]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    else
        throw(ArgumentError("Unknown type_fuel_cell: $type_fuel_cell"))
    end

    return i_exp_cali_t .* 1e4, U_exp_cali_t  # Conversion in A.m-2
end


"""
    plot_experimental_polarisation_curve(type_fuel_cell, i_fc_t, U_exp_t, ax)

This function plots the experimental polarisation curve on the same graph as the model results.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell used in the model. This parameter includes the fuel cell used
  in the model and the corresponding operating conditions.
- `i_fc_t::Vector{<:Real}`: Current density values.
- `U_exp_t::Vector{<:Real}`: Experimental values of the voltage.
- `ax`: Axes-like object exposing `scatter` and `legend` methods.
"""
function plot_experimental_polarisation_curve(type_fuel_cell::String,
                                              i_fc_t::Vector{<:Real},
                                              U_exp_t::Vector{<:Real},
                                              ax)
    # ZSW-GenStack
    if type_fuel_cell == "ZSW-GenStack"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="s", color="black", label="Exp. - nominal")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="v", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 1.61/1.41 bar")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="^", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 2.01/1.81 bar")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="p", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 2.4/2.2 bar")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="D", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 2.8/2.6 bar")
    elseif type_fuel_cell == "ZSW-GenStack_T_62"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="P", color="black", label="Exp. - T = 62 \$^\\circ\$C")
    elseif type_fuel_cell == "ZSW-GenStack_T_76"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="X", color="black", label="Exp. - T = 76 \$^\\circ\$C")
    elseif type_fuel_cell == "ZSW-GenStack_T_84"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="*", color="black", label="Exp. - T = 84 \$^\\circ\$C")

    # EH-31
    elseif type_fuel_cell == "EH-31_1.5"  # at 1.5 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="s", color="black", label="Exp. - P = 1.5 bar")
    elseif type_fuel_cell == "EH-31_2.0"  # at 2.0 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="v", color="black", label="Exp. - P = 2.0 bar")
    elseif type_fuel_cell == "EH-31_2.25"  # at 2.25 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="^", color="black", label="Exp. - P = 2.25 bar")
    elseif type_fuel_cell == "EH-31_2.5"  # at 2.5 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="p", color="black", label="Exp. - P = 2.5 bar")
    else
        throw(ArgumentError("Unknown type_fuel_cell: $type_fuel_cell"))
    end

    ax.legend(loc="best", markerscale=0.5)
    return nothing
end

