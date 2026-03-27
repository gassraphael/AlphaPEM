# -*- coding: utf-8 -*-

"""This file is designated for experimental polarization-curve data used by AlphaPEM."""

# ___________________________________________________Experimental data__________________________________________________

"""
    pola_exp_values(type_fuel_cell, voltage_zone)

This function returns the experimental values of polarisation curves made on different fuel cells at different
operating conditions. The experimental values are used to compare the model results with the experimental data.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell used in the model. This parameter includes the fuel cell used
  in the model and the corresponding operating conditions.
- `voltage_zone::String`: Zone of the polarization curve which is considered. It can be `"full"` or
  `"before_voltage_drop"`.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: Experimental values of current density and voltage.
"""
function pola_exp_values(type_fuel_cell::String,
                         voltage_zone::String)::Tuple{Vector{Float64}, Vector{Float64}}

    # GenStack is a fuel cell developed in open source by ZSW (https://zenodo.org/records/14223364).
    if type_fuel_cell == "ZSW-GenStack"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.200, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp_t = [0.953, 0.864, 0.838, 0.819, 0.804, 0.778, 0.760, 0.743, 0.721, 0.703,
                       0.694, 0.685, 0.666, 0.644, 0.620, 0.593, 0.579, 0.546, 0.486]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.200, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700]
            U_exp_t = [0.953, 0.864, 0.838, 0.819, 0.804, 0.778, 0.760, 0.743, 0.721, 0.703,
                       0.694, 0.685, 0.666, 0.644, 0.620]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700, 0.900,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.171]
            U_exp_t = [0.936, 0.835, 0.809, 0.795, 0.783, 0.759, 0.741, 0.725, 0.701, 0.670,
                       0.661, 0.633, 0.587, 0.541, 0.500, 0.457, 0.437, 0.402]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700]
            U_exp_t = [0.936, 0.835, 0.809, 0.795, 0.783, 0.759, 0.741, 0.725, 0.701]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.415]
            U_exp_t = [0.946, 0.855, 0.830, 0.811, 0.795, 0.770, 0.752, 0.736, 0.717, 0.697,
                       0.685, 0.677, 0.655, 0.629, 0.599, 0.564, 0.545, 0.502, 0.450]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300]
            U_exp_t = [0.946, 0.855, 0.830, 0.811, 0.795, 0.770, 0.752, 0.736, 0.717, 0.697,
                       0.685, 0.677, 0.655]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp_t = [0.949, 0.867, 0.841, 0.821, 0.807, 0.781, 0.763, 0.746, 0.725, 0.706,
                       0.697, 0.687, 0.670, 0.651, 0.630, 0.607, 0.595, 0.566, 0.514]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900]
            U_exp_t = [0.949, 0.867, 0.841, 0.821, 0.807, 0.781, 0.763, 0.746, 0.725, 0.706,
                       0.697, 0.687, 0.670, 0.651, 0.630, 0.607]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp_t = [0.947, 0.872, 0.846, 0.827, 0.812, 0.787, 0.768, 0.752, 0.731, 0.711,
                       0.703, 0.694, 0.676, 0.659, 0.641, 0.622, 0.610, 0.588, 0.547]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900]
            U_exp_t = [0.947, 0.872, 0.846, 0.827, 0.812, 0.787, 0.768, 0.752, 0.731, 0.711,
                       0.703, 0.694, 0.676, 0.659, 0.641, 0.622]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_T_62"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp_t = [0.944, 0.855, 0.827, 0.808, 0.795, 0.771, 0.754, 0.739, 0.717, 0.696,
                       0.685, 0.675, 0.653, 0.631, 0.606, 0.581, 0.566, 0.532, 0.471]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500]
            U_exp_t = [0.944, 0.855, 0.827, 0.808, 0.795, 0.771, 0.754, 0.739, 0.717, 0.696,
                       0.685, 0.675, 0.653, 0.631]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_T_76"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.900,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp_t = [0.946, 0.849, 0.825, 0.811, 0.799, 0.776, 0.759, 0.744, 0.724, 0.702,
                       0.691, 0.679, 0.652, 0.621, 0.587, 0.547, 0.527, 0.482, 0.406]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.900,
                       0.999, 1.099]
            U_exp_t = [0.946, 0.849, 0.825, 0.811, 0.799, 0.776, 0.759, 0.744, 0.724, 0.702,
                       0.691, 0.679]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "ZSW-GenStack_T_84"
        if voltage_zone == "full"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.300, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000]
            U_exp_t = [0.930, 0.847, 0.820, 0.805, 0.794, 0.772, 0.756, 0.741, 0.718, 0.686,
                       0.668, 0.650, 0.614, 0.575, 0.532, 0.486, 0.461]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.001, 0.050, 0.099, 0.150, 0.199, 0.300, 0.400, 0.498, 0.700]
            U_exp_t = [0.930, 0.847, 0.820, 0.805, 0.794, 0.772, 0.756, 0.741, 0.718]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    # EH-31 fuel cell
    elseif type_fuel_cell == "EH-31_1.5"  # at 1.5 bar
        if voltage_zone == "full"
            i_exp_t = [0.050, 0.068, 0.089, 0.110, 0.147, 0.185, 0.233, 0.293, 0.352, 0.395,
                       0.455, 0.510, 0.556, 0.620, 0.672, 0.738, 0.799, 0.850, 0.892, 0.942,
                       1.039, 1.139, 1.212, 1.269, 1.360, 1.432, 1.525, 1.604, 1.683, 1.765,
                       1.878, 1.966, 2.050, 2.109, 2.151, 2.188, 2.246]
            U_exp_t = [0.900, 0.882, 0.865, 0.850, 0.834, 0.823, 0.811, 0.794, 0.781, 0.772,
                       0.761, 0.752, 0.745, 0.735, 0.728, 0.719, 0.712, 0.706, 0.700, 0.694,
                       0.681, 0.668, 0.660, 0.653, 0.641, 0.634, 0.622, 0.610, 0.599, 0.586,
                       0.570, 0.556, 0.540, 0.530, 0.521, 0.513, 0.500]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.050, 0.068, 0.089, 0.110, 0.147, 0.185, 0.233, 0.293, 0.352, 0.395,
                       0.455, 0.510, 0.556, 0.620, 0.672, 0.738, 0.799, 0.850, 0.892, 0.942,
                       1.039, 1.139, 1.212, 1.269, 1.360, 1.432, 1.525, 1.604, 1.683]
            U_exp_t = [0.900, 0.882, 0.865, 0.850, 0.834, 0.823, 0.811, 0.794, 0.781, 0.772,
                       0.761, 0.752, 0.745, 0.735, 0.728, 0.719, 0.712, 0.706, 0.700, 0.694,
                       0.681, 0.668, 0.660, 0.653, 0.641, 0.634, 0.622, 0.610, 0.599]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "EH-31_2.0"  # at 2.0 bar
        if voltage_zone == "full"
            i_exp_t = [0.050, 0.057, 0.079, 0.106, 0.135, 0.171, 0.206, 0.242, 0.302, 0.346,
                       0.395, 0.434, 0.476, 0.531, 0.570, 0.623, 0.681, 0.731, 0.779, 0.822,
                       0.868, 0.930, 0.976, 1.031, 1.090, 1.134, 1.205, 1.242, 1.312, 1.358,
                       1.403, 1.453, 1.501, 1.569, 1.634, 1.725, 1.786, 1.857, 1.924, 1.979,
                       2.050, 2.125, 2.168, 2.214, 2.258, 2.308, 2.348, 2.413, 2.459]
            U_exp_t = [0.900, 0.889, 0.874, 0.860, 0.853, 0.845, 0.837, 0.830, 0.817, 0.808,
                       0.800, 0.792, 0.786, 0.779, 0.772, 0.765, 0.759, 0.753, 0.747, 0.742,
                       0.737, 0.730, 0.726, 0.720, 0.714, 0.710, 0.702, 0.698, 0.690, 0.684,
                       0.679, 0.673, 0.668, 0.659, 0.651, 0.640, 0.631, 0.620, 0.608, 0.598,
                       0.586, 0.573, 0.565, 0.557, 0.548, 0.537, 0.528, 0.513, 0.502]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.050, 0.057, 0.079, 0.106, 0.135, 0.171, 0.206, 0.242, 0.302, 0.346,
                       0.395, 0.434, 0.476, 0.531, 0.570, 0.623, 0.681, 0.731, 0.779, 0.822,
                       0.868, 0.930, 0.976, 1.031, 1.090, 1.134, 1.205, 1.242]
            U_exp_t = [0.900, 0.889, 0.874, 0.860, 0.853, 0.845, 0.837, 0.830, 0.817, 0.808,
                       0.800, 0.792, 0.786, 0.779, 0.772, 0.765, 0.759, 0.753, 0.747, 0.742,
                       0.737, 0.730, 0.726, 0.720, 0.714, 0.710, 0.702, 0.698]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "EH-31_2.25"  # at 2.25 bar
        if voltage_zone == "full"
            i_exp_t = [0.056, 0.095, 0.120, 0.138, 0.160, 0.183, 0.218, 0.248, 0.279, 0.315,
                       0.364, 0.409, 0.477, 0.536, 0.594, 0.641, 0.697, 0.748, 0.809, 0.866,
                       0.944, 1.011, 1.074, 1.142, 1.193, 1.252, 1.322, 1.381, 1.442, 1.496,
                       1.545, 1.599, 1.675, 1.746, 1.827, 1.868, 1.918, 2.004, 2.053, 2.114,
                       2.156, 2.209, 2.257, 2.310, 2.356, 2.403, 2.468, 2.513, 2.552, 2.600,
                       2.636, 2.679, 2.728, 2.794]
            U_exp_t = [0.894, 0.882, 0.873, 0.867, 0.861, 0.854, 0.847, 0.840, 0.834, 0.827,
                       0.819, 0.812, 0.801, 0.793, 0.786, 0.781, 0.775, 0.771, 0.764, 0.759,
                       0.751, 0.746, 0.740, 0.734, 0.728, 0.723, 0.715, 0.709, 0.703, 0.698,
                       0.692, 0.686, 0.678, 0.670, 0.660, 0.654, 0.647, 0.635, 0.628, 0.618,
                       0.613, 0.604, 0.596, 0.587, 0.580, 0.570, 0.559, 0.551, 0.545, 0.536,
                       0.528, 0.520, 0.511, 0.497]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.056, 0.095, 0.120, 0.138, 0.160, 0.183, 0.218, 0.248, 0.279, 0.315,
                       0.364, 0.409, 0.477, 0.536, 0.594, 0.641, 0.697, 0.748, 0.809, 0.866,
                       0.944, 1.011, 1.074, 1.142, 1.193, 1.252, 1.322, 1.381, 1.442, 1.496,
                       1.545, 1.599, 1.675]
            U_exp_t = [0.894, 0.882, 0.873, 0.867, 0.861, 0.854, 0.847, 0.840, 0.834, 0.827,
                       0.819, 0.812, 0.801, 0.793, 0.786, 0.781, 0.775, 0.771, 0.764, 0.759,
                       0.751, 0.746, 0.740, 0.734, 0.728, 0.723, 0.715, 0.709, 0.703, 0.698,
                       0.692, 0.686, 0.678]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    elseif type_fuel_cell == "EH-31_2.5"  # at 2.5 bar
        if voltage_zone == "full"
            i_exp_t = [0.057, 0.070, 0.082, 0.101, 0.127, 0.145, 0.168, 0.200, 0.234, 0.267,
                       0.296, 0.331, 0.355, 0.388, 0.423, 0.467, 0.527, 0.577, 0.632, 0.685,
                       0.740, 0.789, 0.845, 0.898, 0.953, 1.030, 1.124, 1.192, 1.254, 1.314,
                       1.364, 1.434, 1.514, 1.587, 1.643, 1.707, 1.769, 1.826, 1.892, 1.972,
                       2.040, 2.124, 2.192, 2.265, 2.358, 2.429, 2.508, 2.572, 2.624, 2.691,
                       2.750, 2.822, 2.879, 2.918, 2.956, 2.988]
            U_exp_t = [0.900, 0.892, 0.884, 0.875, 0.866, 0.861, 0.856, 0.850, 0.845, 0.840,
                       0.835, 0.829, 0.824, 0.820, 0.814, 0.807, 0.800, 0.793, 0.787, 0.783,
                       0.778, 0.775, 0.771, 0.767, 0.763, 0.758, 0.750, 0.744, 0.738, 0.732,
                       0.726, 0.719, 0.712, 0.703, 0.697, 0.691, 0.685, 0.679, 0.672, 0.663,
                       0.657, 0.648, 0.640, 0.632, 0.621, 0.610, 0.600, 0.591, 0.584, 0.575,
                       0.566, 0.555, 0.546, 0.537, 0.531, 0.524]
        elseif voltage_zone == "before_voltage_drop"
            i_exp_t = [0.057, 0.070, 0.082, 0.101, 0.127, 0.145, 0.168, 0.200, 0.234, 0.267,
                       0.296, 0.331, 0.355, 0.388, 0.423, 0.467, 0.527, 0.577, 0.632, 0.685,
                       0.740, 0.789, 0.845, 0.898, 0.953, 1.030, 1.124, 1.192, 1.254, 1.314,
                       1.364, 1.434, 1.514]
            U_exp_t = [0.900, 0.892, 0.884, 0.875, 0.866, 0.861, 0.856, 0.850, 0.845, 0.840,
                       0.835, 0.829, 0.824, 0.820, 0.814, 0.807, 0.800, 0.793, 0.787, 0.783,
                       0.778, 0.775, 0.771, 0.767, 0.763, 0.758, 0.750, 0.744, 0.738, 0.732,
                       0.726, 0.719, 0.712]
        else
            throw(ArgumentError("The voltage_zone should be either 'full' or 'before_voltage_drop'."))
        end

    else
        throw(ArgumentError("Unknown type_fuel_cell: $type_fuel_cell"))
    end

    return i_exp_t .* 1e4, U_exp_t  # Conversion in A.m-2
end


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

