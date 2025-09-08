# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np


# ___________________________________________________Experimental data__________________________________________________

def pola_exp_values(type_fuel_cell):
    """
    This function returns the experimental values of polarisation curves made on different fuel cells at different
    operating conditions. The experimental values are used to compare the model results with the experimental data.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell used in the model. This parameter includes the fuel cell used in the model and the
        corresponding operating conditions.

    Returns
    -------
    i_exp_t : numpy.ndarray
        Experimental values of the current density.
    U_exp_t : numpy.ndarray
        Experimental values of the voltage.

    """
    # ZSW fuel cell
    if type_fuel_cell == "ZSW-GenStack":
        # # Current density
        # i_exp_t = np.zeros(19)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.200
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18] = 1.900, 2.000, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(19)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.953, 0.864, 0.838, 0.819, 0.804
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.778, 0.760, 0.743, 0.721, 0.703
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.694, 0.685, 0.666, 0.644, 0.620
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18] = 0.593, 0.579, 0.546, 0.486
        # Current density
        i_exp_t = np.zeros(15)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.200
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # Voltage
        U_exp_t = np.zeros(15)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.953, 0.864, 0.838, 0.819, 0.804
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.778, 0.760, 0.743, 0.721, 0.703
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.694, 0.685, 0.666, 0.644, 0.620

    elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41":
        # # Current density
        # i_exp_t = np.zeros(18)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.200
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.300, 0.400, 0.498, 0.700, 0.900
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16], i_exp_t[17] = 1.900, 2.000, 2.171
        # # Voltage
        # U_exp_t = np.zeros(18)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.936, 0.835, 0.809, 0.795, 0.783
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.759, 0.741, 0.725, 0.701, 0.670
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.661, 0.633, 0.587, 0.541, 0.500
        # U_exp_t[15], U_exp_t[16], U_exp_t[17] = 0.457, 0.437, 0.402
        # Current density
        i_exp_t = np.zeros(9)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.200
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8] = 0.300, 0.400, 0.498, 0.700
        # Voltage
        U_exp_t = np.zeros(9)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.936, 0.835, 0.809, 0.795, 0.783
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8] = 0.759, 0.741, 0.725, 0.701

    elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81":
        # # Current density
        # i_exp_t = np.zeros(19)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18] = 1.900, 2.000, 2.200, 2.415
        # # Voltage
        # U_exp_t = np.zeros(19)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.946, 0.855, 0.830, 0.811, 0.795
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.770, 0.752, 0.736, 0.717, 0.697
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.685, 0.677, 0.655, 0.629, 0.599
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18] = 0.564, 0.545, 0.502, 0.450
        # Current density
        i_exp_t = np.zeros(13)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        i_exp_t[10], i_exp_t[11], i_exp_t[12] = 0.999, 1.099, 1.300
        # Voltage
        U_exp_t = np.zeros(13)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.946, 0.855, 0.830, 0.811, 0.795
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.770, 0.752, 0.736, 0.717, 0.697
        U_exp_t[10], U_exp_t[11], U_exp_t[12] = 0.685, 0.677, 0.655

    elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2":
        # # Current density
        # i_exp_t = np.zeros(19)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.200
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.300, 0.400, 0.498, 0.700, 0.901
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18] = 1.900, 2.000, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(19)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.949, 0.867, 0.841, 0.821, 0.807
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.781, 0.763, 0.746, 0.725, 0.706
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.697, 0.687, 0.670, 0.651, 0.630
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18] = 0.607, 0.595, 0.566, 0.514
        # Current density
        i_exp_t = np.zeros(16)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.200
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.300, 0.400, 0.498, 0.700, 0.901
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        i_exp_t[15] = 1.900
        # Voltage
        U_exp_t = np.zeros(16)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.949, 0.867, 0.841, 0.821, 0.807
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.781, 0.763, 0.746, 0.725, 0.706
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.697, 0.687, 0.670, 0.651, 0.630
        U_exp_t[15] = 0.607

    elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6":
        # # Current density
        # i_exp_t = np.zeros(19)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18] = 1.900, 2.000, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(19)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.947, 0.872, 0.846, 0.827, 0.812
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.787, 0.768, 0.752, 0.731, 0.711
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.703, 0.694, 0.676, 0.659, 0.641
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18] = 0.622, 0.610, 0.588, 0.547
        # Current density
        i_exp_t = np.zeros(16)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        i_exp_t[15] = 1.900
        # Voltage
        U_exp_t = np.zeros(16)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.947, 0.872, 0.846, 0.827, 0.812
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.787, 0.768, 0.752, 0.731, 0.711
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.703, 0.694, 0.676, 0.659, 0.641
        U_exp_t[15] = 0.622

    elif type_fuel_cell == "ZSW-GenStack_T_62":
        # # Current density
        # i_exp_t = np.zeros(19)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18] = 1.900, 2.000, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(19)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.944, 0.855, 0.827, 0.808, 0.795
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.771, 0.754, 0.739, 0.717, 0.696
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.685, 0.675, 0.653, 0.631, 0.606
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18] = 0.581, 0.566, 0.532, 0.471
        # Current density
        i_exp_t = np.zeros(14)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.901
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13] = 0.999, 1.099, 1.300, 1.500
        # Voltage
        U_exp_t = np.zeros(14)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.944, 0.855, 0.827, 0.808, 0.795
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.771, 0.754, 0.739, 0.717, 0.696
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13] = 0.685, 0.675, 0.653, 0.631

    elif type_fuel_cell == "ZSW-GenStack_T_76":
        # # Current density
        # i_exp_t = np.zeros(19)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.900
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18] = 1.900, 2.000, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(19)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.946, 0.849, 0.825, 0.811, 0.799
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.776, 0.759, 0.744, 0.724, 0.702
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.691, 0.679, 0.652, 0.621, 0.587
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18] = 0.547, 0.527, 0.482, 0.406
        # Current density
        i_exp_t = np.zeros(12)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.299, 0.400, 0.498, 0.700, 0.900
        i_exp_t[10], i_exp_t[11] = 0.999, 1.099
        # Voltage
        U_exp_t = np.zeros(12)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.946, 0.849, 0.825, 0.811, 0.799
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.776, 0.759, 0.744, 0.724, 0.702
        U_exp_t[10], U_exp_t[11] = 0.691, 0.679

    elif type_fuel_cell == "ZSW-GenStack_T_84":
        # # Current density
        # i_exp_t = np.zeros(17)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.300, 0.400, 0.498, 0.700, 0.901
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.999, 1.099, 1.300, 1.500, 1.700
        # i_exp_t[15], i_exp_t[16] = 1.900, 2.000
        # # Voltage
        # U_exp_t = np.zeros(17)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.930, 0.847, 0.820, 0.805, 0.794
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.772, 0.756, 0.741, 0.718, 0.686
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.668, 0.650, 0.614, 0.575, 0.532
        # U_exp_t[15], U_exp_t[16] = 0.486, 0.461
        # Current density
        i_exp_t = np.zeros(9)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.099, 0.150, 0.199
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8] = 0.300, 0.400, 0.498, 0.700
        # Voltage
        U_exp_t = np.zeros(9)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.930, 0.847, 0.820, 0.805, 0.794
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8] = 0.772, 0.756, 0.741, 0.718

    # EH-31 fuel cell
    elif type_fuel_cell == "EH-31_1.5":  # at 1.5 bar
        # Current density
        i_exp_t = np.zeros(37)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.050, 0.068, 0.089, 0.110, 0.147
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.185, 0.233, 0.293, 0.352, 0.395
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.455, 0.510, 0.556, 0.620, 0.672
        i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.738, 0.799, 0.850, 0.892, 0.942
        i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 1.039, 1.139, 1.212, 1.269, 1.360
        i_exp_t[25], i_exp_t[26], i_exp_t[27], i_exp_t[28], i_exp_t[29] = 1.432, 1.525, 1.604, 1.683, 1.765
        i_exp_t[30], i_exp_t[31], i_exp_t[32], i_exp_t[33], i_exp_t[34] = 1.878, 1.966, 2.050, 2.109, 2.151
        i_exp_t[35], i_exp_t[36] = 2.188, 2.246
        # Voltage
        U_exp_t = np.zeros(37)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.900, 0.882, 0.865, 0.850, 0.834
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.823, 0.811, 0.794, 0.781, 0.772
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.761, 0.752, 0.745, 0.735, 0.728
        U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.719, 0.712, 0.706, 0.700, 0.694
        U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.681, 0.668, 0.660, 0.653, 0.641
        U_exp_t[25], U_exp_t[26], U_exp_t[27], U_exp_t[28], U_exp_t[29] = 0.634, 0.622, 0.610, 0.599, 0.586
        U_exp_t[30], U_exp_t[31], U_exp_t[32], U_exp_t[33], U_exp_t[34] = 0.570, 0.556, 0.540, 0.530, 0.521
        U_exp_t[35], U_exp_t[36] = 0.513, 0.500
        # # Current density
        # i_exp_t = np.zeros(29)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.050, 0.068, 0.089, 0.110, 0.147
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.185, 0.233, 0.293, 0.352, 0.395
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.455, 0.510, 0.556, 0.620, 0.672
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.738, 0.799, 0.850, 0.892, 0.942
        # i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 1.039, 1.139, 1.212, 1.269, 1.360
        # i_exp_t[25], i_exp_t[26], i_exp_t[27], i_exp_t[28] = 1.432, 1.525, 1.604, 1.683
        # # Voltage
        # U_exp_t = np.zeros(29)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.900, 0.882, 0.865, 0.850, 0.834
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.823, 0.811, 0.794, 0.781, 0.772
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.761, 0.752, 0.745, 0.735, 0.728
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.719, 0.712, 0.706, 0.700, 0.694
        # U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.681, 0.668, 0.660, 0.653, 0.641
        # U_exp_t[25], U_exp_t[26], U_exp_t[27], U_exp_t[28] = 0.634, 0.622, 0.610, 0.599

    elif type_fuel_cell == "EH-31_2.0":  # at 2.0 bar
        # Current density
        i_exp_t = np.zeros(49)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.050, 0.057, 0.079, 0.106, 0.135
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.171, 0.206, 0.242, 0.302, 0.346
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.395, 0.434, 0.476, 0.531, 0.570
        i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.623, 0.681, 0.731, 0.779, 0.822
        i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 0.868, 0.930, 0.976, 1.031, 1.090
        i_exp_t[25], i_exp_t[26], i_exp_t[27], i_exp_t[28], i_exp_t[29] = 1.134, 1.205, 1.242, 1.312, 1.358
        i_exp_t[30], i_exp_t[31], i_exp_t[32], i_exp_t[33], i_exp_t[34] = 1.403, 1.453, 1.501, 1.569, 1.634,
        i_exp_t[35], i_exp_t[36], i_exp_t[37], i_exp_t[38], i_exp_t[39] = 1.725, 1.786, 1.857, 1.924, 1.979
        i_exp_t[40], i_exp_t[41], i_exp_t[42], i_exp_t[43], i_exp_t[44] = 2.050, 2.125, 2.168, 2.214, 2.258
        i_exp_t[45], i_exp_t[46], i_exp_t[47], i_exp_t[48] = 2.308, 2.348, 2.413, 2.459
        # Voltage
        U_exp_t = np.zeros(49)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.900, 0.889, 0.874, 0.860, 0.853
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.845, 0.837, 0.830, 0.817, 0.808
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.800, 0.792, 0.786, 0.779, 0.772
        U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.765, 0.759, 0.753, 0.747, 0.742,
        U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.737, 0.730, 0.726, 0.720, 0.714,
        U_exp_t[25], U_exp_t[26], U_exp_t[27], U_exp_t[28], U_exp_t[29] = 0.710, 0.702, 0.698, 0.690, 0.684
        U_exp_t[30], U_exp_t[31], U_exp_t[32], U_exp_t[33], U_exp_t[34] = 0.679, 0.673, 0.668, 0.659, 0.651
        U_exp_t[35], U_exp_t[36], U_exp_t[37], U_exp_t[38], U_exp_t[39] = 0.640, 0.631, 0.620, 0.608, 0.598
        U_exp_t[40], U_exp_t[41], U_exp_t[42], U_exp_t[43], U_exp_t[44] = 0.586, 0.573, 0.565, 0.557, 0.548
        U_exp_t[45], U_exp_t[46], U_exp_t[47], U_exp_t[48] = 0.537, 0.528, 0.513, 0.502
        # # Current density
        # i_exp_t = np.zeros(28)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.050, 0.057, 0.079, 0.106, 0.135
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.171, 0.206, 0.242, 0.302, 0.346
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.395, 0.434, 0.476, 0.531, 0.570
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.623, 0.681, 0.731, 0.779, 0.822
        # i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 0.868, 0.930, 0.976, 1.031, 1.090
        # i_exp_t[25], i_exp_t[26], i_exp_t[27] = 1.134, 1.205, 1.242
        # # Voltage
        # U_exp_t = np.zeros(28)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.900, 0.889, 0.874, 0.860, 0.853
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.845, 0.837, 0.830, 0.817, 0.808
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.800, 0.792, 0.786, 0.779, 0.772
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.765, 0.759, 0.753, 0.747, 0.742,
        # U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.737, 0.730, 0.726, 0.720, 0.714,
        # U_exp_t[25], U_exp_t[26], U_exp_t[27] = 0.710, 0.702, 0.698

    elif type_fuel_cell == "EH-31_2.25":  # at 2.25 bar
        # Current density
        i_exp_t = np.zeros(54)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.056, 0.095, 0.120, 0.138, 0.160
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.183, 0.218, 0.248, 0.279, 0.315
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.364, 0.409, 0.477, 0.536, 0.594
        i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.641, 0.697, 0.748, 0.809, 0.866
        i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 0.944, 1.011, 1.074, 1.142, 1.193
        i_exp_t[25], i_exp_t[26], i_exp_t[27], i_exp_t[28], i_exp_t[29] = 1.252, 1.322, 1.381, 1.442, 1.496
        i_exp_t[30], i_exp_t[31], i_exp_t[32], i_exp_t[33], i_exp_t[34] = 1.545, 1.599, 1.675, 1.746, 1.827
        i_exp_t[35], i_exp_t[36], i_exp_t[37], i_exp_t[38], i_exp_t[39] = 1.868, 1.918, 2.004, 2.053, 2.114
        i_exp_t[40], i_exp_t[41], i_exp_t[42], i_exp_t[43], i_exp_t[44] = 2.156, 2.209, 2.257, 2.310, 2.356
        i_exp_t[45], i_exp_t[46], i_exp_t[47], i_exp_t[48], i_exp_t[49] = 2.403, 2.468, 2.513, 2.552, 2.600
        i_exp_t[50], i_exp_t[51], i_exp_t[52], i_exp_t[53] = 2.636, 2.679, 2.728, 2.794
        # Voltage
        U_exp_t = np.zeros(54)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.894, 0.882, 0.873, 0.867, 0.861
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.854, 0.847, 0.840, 0.834, 0.827
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.819, 0.812, 0.801, 0.793, 0.786
        U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.781, 0.775, 0.771, 0.764, 0.759
        U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.751, 0.746, 0.740, 0.734, 0.728
        U_exp_t[25], U_exp_t[26], U_exp_t[27], U_exp_t[28], U_exp_t[29] = 0.723, 0.715, 0.709, 0.703, 0.698
        U_exp_t[30], U_exp_t[31], U_exp_t[32], U_exp_t[33], U_exp_t[34] = 0.692, 0.686, 0.678, 0.670, 0.660
        U_exp_t[35], U_exp_t[36], U_exp_t[37], U_exp_t[38], U_exp_t[39] = 0.654, 0.647, 0.635, 0.628, 0.618
        U_exp_t[40], U_exp_t[41], U_exp_t[42], U_exp_t[43], U_exp_t[44] = 0.613, 0.604, 0.596, 0.587, 0.580
        U_exp_t[45], U_exp_t[46], U_exp_t[47], U_exp_t[48], U_exp_t[49] = 0.570, 0.559, 0.551, 0.545, 0.536
        U_exp_t[50], U_exp_t[51], U_exp_t[52], U_exp_t[53] = 0.528, 0.520, 0.511, 0.497
        # # Current density
        # i_exp_t = np.zeros(33)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.056, 0.095, 0.120, 0.138, 0.160
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.183, 0.218, 0.248, 0.279, 0.315
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.364, 0.409, 0.477, 0.536, 0.594
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.641, 0.697, 0.748, 0.809, 0.866
        # i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 0.944, 1.011, 1.074, 1.142, 1.193
        # i_exp_t[25], i_exp_t[26], i_exp_t[27], i_exp_t[28], i_exp_t[29] = 1.252, 1.322, 1.381, 1.442, 1.496
        # i_exp_t[30], i_exp_t[31], i_exp_t[32] = 1.545, 1.599, 1.675
        # # Voltage
        # U_exp_t = np.zeros(33)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.894, 0.882, 0.873, 0.867, 0.861
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.854, 0.847, 0.840, 0.834, 0.827
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.819, 0.812, 0.801, 0.793, 0.786
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.781, 0.775, 0.771, 0.764, 0.759
        # U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.751, 0.746, 0.740, 0.734, 0.728
        # U_exp_t[25], U_exp_t[26], U_exp_t[27], U_exp_t[28], U_exp_t[29] = 0.723, 0.715, 0.709, 0.703, 0.698
        # U_exp_t[30], U_exp_t[31], U_exp_t[32] = 0.692, 0.686, 0.678

    elif type_fuel_cell == "EH-31_2.5":  # at 2.5 bar
        # Current density
        i_exp_t = np.zeros(56)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.057, 0.070, 0.082, 0.101, 0.127
        i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.145, 0.168, 0.200, 0.234, 0.267
        i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.296, 0.331, 0.355, 0.388, 0.423
        i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.467, 0.527, 0.577, 0.632, 0.685
        i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 0.740, 0.789, 0.845, 0.898, 0.953
        i_exp_t[25], i_exp_t[26], i_exp_t[27], i_exp_t[28], i_exp_t[29] = 1.030, 1.124, 1.192, 1.254, 1.314
        i_exp_t[30], i_exp_t[31], i_exp_t[32], i_exp_t[33], i_exp_t[34] = 1.364, 1.434, 1.514, 1.587, 1.643
        i_exp_t[35], i_exp_t[36], i_exp_t[37], i_exp_t[38], i_exp_t[39] = 1.707, 1.769, 1.826, 1.892, 1.972
        i_exp_t[40], i_exp_t[41], i_exp_t[42], i_exp_t[43], i_exp_t[44] = 2.040, 2.124, 2.192, 2.265, 2.358
        i_exp_t[45], i_exp_t[46], i_exp_t[47], i_exp_t[48], i_exp_t[49] = 2.429, 2.508, 2.572, 2.624, 2.691
        i_exp_t[50], i_exp_t[51], i_exp_t[52], i_exp_t[53], i_exp_t[54] = 2.750, 2.822, 2.879, 2.918, 2.956
        i_exp_t[55] = 2.988
        # Voltage
        U_exp_t = np.zeros(56)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.900, 0.892, 0.884, 0.875, 0.866
        U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.861, 0.856, 0.850, 0.845, 0.840
        U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.835, 0.829, 0.824, 0.820, 0.814
        U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.807, 0.800, 0.793, 0.787, 0.783
        U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.778, 0.775, 0.771, 0.767, 0.763
        U_exp_t[25], U_exp_t[26], U_exp_t[27], U_exp_t[28], U_exp_t[29] = 0.758, 0.750, 0.744, 0.738, 0.732
        U_exp_t[30], U_exp_t[31], U_exp_t[32], U_exp_t[33], U_exp_t[34] = 0.726, 0.719, 0.712, 0.703, 0.697
        U_exp_t[35], U_exp_t[36], U_exp_t[37], U_exp_t[38], U_exp_t[39] = 0.691, 0.685, 0.679, 0.672, 0.663
        U_exp_t[40], U_exp_t[41], U_exp_t[42], U_exp_t[43], U_exp_t[44] = 0.657, 0.648, 0.640, 0.632, 0.621
        U_exp_t[45], U_exp_t[46], U_exp_t[47], U_exp_t[48], U_exp_t[49] = 0.610, 0.600, 0.591, 0.584, 0.575
        U_exp_t[50], U_exp_t[51], U_exp_t[52], U_exp_t[53], U_exp_t[54] = 0.566, 0.555, 0.546, 0.537, 0.531
        U_exp_t[55] = 0.524
        # # Current density
        # i_exp_t = np.zeros(33)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.057, 0.070, 0.082, 0.101, 0.127
        # i_exp_t[5], i_exp_t[6], i_exp_t[7], i_exp_t[8], i_exp_t[9] = 0.145, 0.168, 0.200, 0.234, 0.267
        # i_exp_t[10], i_exp_t[11], i_exp_t[12], i_exp_t[13], i_exp_t[14] = 0.296, 0.331, 0.355, 0.388, 0.423
        # i_exp_t[15], i_exp_t[16], i_exp_t[17], i_exp_t[18], i_exp_t[19] = 0.467, 0.527, 0.577, 0.632, 0.685
        # i_exp_t[20], i_exp_t[21], i_exp_t[22], i_exp_t[23], i_exp_t[24] = 0.740, 0.789, 0.845, 0.898, 0.953
        # i_exp_t[25], i_exp_t[26], i_exp_t[27], i_exp_t[28], i_exp_t[29] = 1.030, 1.124, 1.192, 1.254, 1.314
        # i_exp_t[30], i_exp_t[31], i_exp_t[32] = 1.364, 1.434, 1.514
        # # Voltage
        # U_exp_t = np.zeros(33)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.900, 0.892, 0.884, 0.875, 0.866
        # U_exp_t[5], U_exp_t[6], U_exp_t[7], U_exp_t[8], U_exp_t[9] = 0.861, 0.856, 0.850, 0.845, 0.840
        # U_exp_t[10], U_exp_t[11], U_exp_t[12], U_exp_t[13], U_exp_t[14] = 0.835, 0.829, 0.824, 0.820, 0.814
        # U_exp_t[15], U_exp_t[16], U_exp_t[17], U_exp_t[18], U_exp_t[19] = 0.807, 0.800, 0.793, 0.787, 0.783
        # U_exp_t[20], U_exp_t[21], U_exp_t[22], U_exp_t[23], U_exp_t[24] = 0.778, 0.775, 0.771, 0.767, 0.763
        # U_exp_t[25], U_exp_t[26], U_exp_t[27], U_exp_t[28], U_exp_t[29] = 0.758, 0.750, 0.744, 0.738, 0.732
        # U_exp_t[30], U_exp_t[31], U_exp_t[32] = 0.726, 0.719, 0.712

    return i_exp_t * 1e4, U_exp_t # Conversion in A.m-2

def pola_exp_values_calibration(type_fuel_cell):
    """
    This function returns the experimental values of polarisation curves made on different fuel cells at different
    operating conditions. The experimental values are used for calibrating the model and so are composed of a reduced
    number of points compare to the pola_exp_values function. These points are specifically chosen to be as few as
    possible while still providing a good representation of the polarisation curve.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell used in the model. This parameter includes the fuel cell used in the model and the
        corresponding operating conditions.

    Returns
    -------
    i_exp_t : numpy.ndarray
        Experimental values of the current density.
    U_exp_t : numpy.ndarray
        Experimental values of the voltage.

    """
    # ZSW fuel cell
    if type_fuel_cell == "ZSW-GenStack":
        # # Current density
        # i_exp_t = np.zeros(7)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.498, 1.099
        # i_exp_t[4], i_exp_t[5], i_exp_t[6] = 1.700, 2.000, 2.500
        # # Voltage
        # U_exp_t = np.zeros(7)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.953, 0.864, 0.743, 0.685
        # U_exp_t[4], U_exp_t[5], U_exp_t[6] = 0.620, 0.579, 0.486
        # Current density
        i_exp_t = np.zeros(5)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.498, 1.099, 1.700
        # Voltage
        U_exp_t = np.zeros(5)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.953, 0.864, 0.743, 0.685, 0.620

    elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41":
        # # Current density
        # i_exp_t = np.zeros(7)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.300, 0.700
        # i_exp_t[4], i_exp_t[5], i_exp_t[6] = 0.900, 1.500, 2.171
        # # Voltage
        # U_exp_t = np.zeros(7)
        # U_exp_t[0], U_exp_t[1], U_exp_t[5], U_exp_t[8] = 0.936, 0.835, 0.759, 0.701
        # U_exp_t[9], U_exp_t[13], U_exp_t[17] = 0.670, 0.541, 0.402
        # Current density
        i_exp_t = np.zeros(4)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.300, 0.700
        # Voltage
        U_exp_t = np.zeros(4)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.936, 0.835, 0.759, 0.701

    elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81":
        # # Current density
        # i_exp_t = np.zeros(6)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2] = 0.001, 0.050, 0.498
        # i_exp_t[3], i_exp_t[4], i_exp_t[5] = 1.300, 2.000, 2.415
        # # Voltage
        # U_exp_t = np.zeros(6)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2] = 0.946, 0.855, 0.736
        # U_exp_t[3], U_exp_t[4], U_exp_t[5] = 0.655, 0.545, 0.450
        # Current density
        i_exp_t = np.zeros(4)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.498, 1.300
        # Voltage
        U_exp_t = np.zeros(4)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.946, 0.855, 0.736, 0.655

    elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2":
        # # Current density
        # i_exp_t = np.zeros(7)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.498, 1.099
        # i_exp_t[4], i_exp_t[5], i_exp_t[6] = 1.900, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(7)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.949, 0.867, 0.746, 0.687
        # U_exp_t[4], U_exp_t[5], U_exp_t[6] = 0.607, 0.566, 0.514
        # Current density
        i_exp_t = np.zeros(5)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.498, 1.099, 1.900
        # Voltage
        U_exp_t = np.zeros(5)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.949, 0.867, 0.746, 0.687, 0.607

    elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6":
        # # Current density
        # i_exp_t = np.zeros(7)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.498, 1.099
        # i_exp_t[4], i_exp_t[5], i_exp_t[6] = 1.900, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(7)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.947, 0.872, 0.752, 0.694
        # U_exp_t[4], U_exp_t[5], U_exp_t[6] = 0.622, 0.588, 0.547
        # Current density
        i_exp_t = np.zeros(5)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3], i_exp_t[4] = 0.001, 0.050, 0.498, 1.099, 1.900
        # Voltage
        U_exp_t = np.zeros(5)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3], U_exp_t[4] = 0.947, 0.872, 0.752, 0.703, 0.622

    elif type_fuel_cell == "ZSW-GenStack_T_62":
        # # Current density
        # i_exp_t = np.zeros(6)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2] = 0.001, 0.050, 0.498
        # i_exp_t[3], i_exp_t[4], i_exp_t[5] = 1.500, 2.000, 2.500
        # # Voltage
        # U_exp_t = np.zeros(6)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2] = 0.944, 0.855, 0.739
        # U_exp_t[3], U_exp_t[4], U_exp_t[5] = 0.631, 0.566, 0.471
        # Current density
        i_exp_t = np.zeros(4)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.498, 1.500
        # Voltage
        U_exp_t = np.zeros(4)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.944, 0.855, 0.739, 0.631

    elif type_fuel_cell == "ZSW-GenStack_T_76":
        # # Current density
        # i_exp_t = np.zeros(7)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.498, 1.099
        # i_exp_t[4], i_exp_t[5], i_exp_t[6] = 1.700, 2.200, 2.500
        # # Voltage
        # U_exp_t = np.zeros(7)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.946, 0.849, 0.744, 0.679
        # U_exp_t[4], U_exp_t[5], U_exp_t[6] = 0.587, 0.482, 0.406
        # Current density
        i_exp_t = np.zeros(4)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.498, 1.099
        # Voltage
        U_exp_t = np.zeros(4)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.946, 0.849, 0.744, 0.679

    elif type_fuel_cell == "ZSW-GenStack_T_84":
        # # Current density
        # i_exp_t = np.zeros(7)
        # i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.300, 0.700
        # i_exp_t[4], i_exp_t[5], i_exp_t[6] = 0.901, 1.500, 2.000
        # # Voltage
        # U_exp_t = np.zeros(7)
        # U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.930, 0.847, 0.772, 0.718
        # U_exp_t[4], U_exp_t[5], U_exp_t[6] = 0.686, 0.575, 0.461
        # Current density
        i_exp_t = np.zeros(4)
        i_exp_t[0], i_exp_t[1], i_exp_t[2], i_exp_t[3] = 0.001, 0.050, 0.300, 0.700
        # Voltage
        U_exp_t = np.zeros(4)
        U_exp_t[0], U_exp_t[1], U_exp_t[2], U_exp_t[3] = 0.930, 0.847, 0.772, 0.718

    # EH-31 fuel cell
    elif type_fuel_cell == "EH-31_1.5":  # at 1.5 bar
        # # Current density
        # i_exp_cali_t = np.zeros(7)
        # i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.050, 0.110, 0.293, 1.039
        # i_exp_cali_t[4], i_exp_cali_t[5], i_exp_cali_t[6] = 1.683, 1.966, 2.246
        # # Voltage
        # U_exp_cali_t = np.zeros(7)
        # U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.900, 0.850, 0.794, 0.681
        # U_exp_cali_t[4], U_exp_cali_t[5], U_exp_cali_t[6] = 0.599, 0.556, 0.500
        # Current density
        i_exp_cali_t = np.zeros(5)
        i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.050, 0.110, 0.293, 1.039
        i_exp_cali_t[4] = 1.683
        # Voltage
        U_exp_cali_t = np.zeros(5)
        U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.900, 0.850, 0.794, 0.681
        U_exp_cali_t[4] = 0.599

    elif type_fuel_cell == "EH-31_2.0":  # at 2.0 bar
        # # Current density
        # i_exp_cali_t = np.zeros(8)
        # i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.050, 0.106, 0.242, 0.681
        # i_exp_cali_t[4], i_exp_cali_t[5], i_exp_cali_t[6], i_exp_cali_t[7] = 1.242, 1.501, 1.979, 2.459
        # # Voltage
        # U_exp_cali_t = np.zeros(8)
        # U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.900, 0.860, 0.830, 0.759
        # U_exp_cali_t[4], U_exp_cali_t[5], U_exp_cali_t[6], U_exp_cali_t[7] = 0.698, 0.668, 0.598, 0.502
        # Current density
        i_exp_cali_t = np.zeros(5)
        i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.050, 0.106, 0.242, 0.681
        i_exp_cali_t[4] = 1.242
        # Voltage
        U_exp_cali_t = np.zeros(5)
        U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.900, 0.860, 0.830, 0.759
        U_exp_cali_t[4] = 0.698

    elif type_fuel_cell == "EH-31_2.25":  # at 2.25 bar
        # # Current density
        # i_exp_cali_t = np.zeros(8)
        # i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.056, 0.183, 0.364, 1.011
        # i_exp_cali_t[4], i_exp_cali_t[5], i_exp_cali_t[6], i_exp_cali_t[7] = 1.675, 1.918, 2.356, 2.794
        # # Voltage
        # U_exp_cali_t = np.zeros(8)
        # U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.894, 0.854, 0.819, 0.746
        # U_exp_cali_t[4], U_exp_cali_t[5], U_exp_cali_t[6], U_exp_cali_t[7] = 0.678, 0.647, 0.580, 0.497
        # Current density
        i_exp_cali_t = np.zeros(5)
        i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.056, 0.183, 0.364, 1.011
        i_exp_cali_t[4] = 1.675
        # Voltage
        U_exp_cali_t = np.zeros(5)
        U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.894, 0.854, 0.819, 0.746
        U_exp_cali_t[4] = 0.678

    elif type_fuel_cell == "EH-31_2.5":  # at 2.5 bar
        # # Current density
        # i_exp_cali_t = np.zeros(10)
        # i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.057, 0.127, 0.296, 0.527
        # i_exp_cali_t[4], i_exp_cali_t[5], i_exp_cali_t[6], i_exp_cali_t[7] = 1.030, 1.514, 1.972, 2.358
        # i_exp_cali_t[8], i_exp_cali_t[9] = 2.691, 2.988
        # # Voltage
        # U_exp_cali_t = np.zeros(10)
        # U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.900, 0.866, 0.835, 0.800
        # U_exp_cali_t[4], U_exp_cali_t[5], U_exp_cali_t[6], U_exp_cali_t[7] = 0.758, 0.712, 0.663, 0.621
        # U_exp_cali_t[8], U_exp_cali_t[9] = 0.575, 0.524
        # Current density
        i_exp_cali_t = np.zeros(6)
        i_exp_cali_t[0], i_exp_cali_t[1], i_exp_cali_t[2], i_exp_cali_t[3] = 0.057, 0.127, 0.296, 0.527
        i_exp_cali_t[4], i_exp_cali_t[5] = 1.030, 1.514
        # Voltage
        U_exp_cali_t = np.zeros(6)
        U_exp_cali_t[0], U_exp_cali_t[1], U_exp_cali_t[2], U_exp_cali_t[3] = 0.900, 0.866, 0.835, 0.800
        U_exp_cali_t[4], U_exp_cali_t[5] = 0.758, 0.712

    return i_exp_cali_t * 1e4, U_exp_cali_t # Conversion in A.m-2

def plot_experimental_polarisation_curve(type_fuel_cell, i_fc_t, U_exp_t, ax):
    """
    This function plots the experimental polarisation curve on the same graph as the model results.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell used in the model. This parameter includes the fuel cell used in the model and the
        corresponding operating conditions.
    i_fc_t : numpy.ndarray
        Current density values.
    U_exp_t : numpy.ndarray
        Experimental values of the voltage.
    ax : matplotlib.axes.Axes
        Axes object on which the experimental data is plotted.
    """
    # ZSW-GenStack
    if type_fuel_cell == "ZSW-GenStack":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - nominal operating conditions")
    elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - P$_a$ = 1.61 bar, P$_c$ = 1.41 bar")
    elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - P$_a$ = 2.01 bar, P$_c$ = 1.81 bar")
    elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - P$_a$ = 2.4 bar, P$_c$ = 2.2 bar")
    elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - P$_a$ = 2.8 bar, P$_c$ = 2.6 bar")
    elif type_fuel_cell == "ZSW-GenStack_T_62":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - T = 62 °C")
    elif type_fuel_cell == "ZSW-GenStack_T_76":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - T = 76 °C")
    elif type_fuel_cell == "ZSW-GenStack_T_84":
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - T = 84 °C")
    # EH-31
    elif type_fuel_cell == "EH-31_1.5":  # at 1.5 bar
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="s", color="black", label="Exp. - P = 1.5 bar")
    elif type_fuel_cell == "EH-31_2.0":  # at 2.0 bar
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="v", color="black", label="Exp. - P = 2.0 bar")
    elif type_fuel_cell == "EH-31_2.25":  # at 2.25 bar
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="^", color="black", label="Exp. - P = 2.25 bar")
    elif type_fuel_cell == "EH-31_2.5":  # at 2.5 bar
        ax.scatter(i_fc_t, U_exp_t, linewidths=1.5, marker="p", color="black", label="Exp. - P = 2.5 bar")

    ax.legend(loc='best', markerscale=0.5)
