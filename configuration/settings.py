# -*- coding: utf-8 -*-

"""This file is used to set the parameters of the fuel cell system.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing functions
from configuration.current_densities import (step_current, polarization_current, polarization_current_for_calibration,
                                             EIS_current)
from modules.settings_modules import stored_operating_inputs, stored_physical_parameters, EIS_parameters

# _______________________________________________________Settings_______________________________________________________
def current_density_parameters(type_current=None):
    """This function is used to set the parameters of the current density which is imposed to the fuel cell system.

    Parameters
    ----------
    type_current : str
        Type of current density which is imposed to the fuel cell system. It can be "step", "polarization" or "EIS".

    Returns
    -------
    step_current_parameters : dict
        Parameters for the step current density. It is a dictionary containing:
        - 'delta_t_ini_step': the initial time (in seconds) at zero current density for the stabilisation of the
        internal states,
        - 'delta_t_load_step': the loading time (in seconds) for the step current density function, from 0 to
        i_step,
        - 'delta_t_break_step': the time (in seconds) at i_step current density for the stabilisation of the
        internal states,
        - 'i_step': the current density (in A.m-2) for the step current density function,
        - 'delta_t_dyn_step': the time (in seconds) for dynamic display of the step current density function.
    pola_current_parameters : dict
        Parameters for the polarization current density. It is a dictionary containing:
        - 'delta_t_ini_pola': the initial time (in seconds) at zero current density for the stabilisation of the
        internal states,
        - 'delta_t_load_pola': the loading time (in seconds) for one step current of the polarisation current
        density function,
        - 'delta_t_break_pola': the breaking time (in seconds) for one step current, for the stabilisation of the
        internal states,
        - 'delta_i_pola': the current density step (in A.m-2) for the polarisation current density function.
        - 'i_max_pola': the maximum current density (in A.m-2) for the polarization curve.
    pola_current_for_cali_parameters : dict
        Parameters for the polarization current density for calibration. It is a dictionary containing:
        - 'delta_t_ini_pola_cali': the initial time (in seconds) at zero current density for the stabilisation of
        the internal states,
        - 'delta_t_load_pola_cali': the loading time (in seconds) for one step current of the polarisation current
        density function,
        - 'delta_t_break_pola_cali': the breaking time (in seconds) for one step current, for the stabilisation of
        the internal states.
    i_EIS : float
        Parameters for the EIS curve. It is the current for which a perturbation is added.
    ratio_EIS : float
        Parameters for the EIS curve. It is the ratio of the current for which a perturbation is added.
    f_EIS : tuple
        Frequency parameters for the EIS_current density function.
        It is a tuple containing the power of the initial frequency 'f_power_min_EIS'
        (f_min_EIS = 10**f_power_min_EIS), the power of the final frequency 'f_power_max_EIS', the number of
        frequencies tested 'nb_f_EIS', and the number of points calculated per specific period 'nb_points_EIS'.
    t_EIS : tuple
        Time parameters for the EIS_current density function.
        It is a tuple containing the initial EIS time after stack equilibrium 't0_EIS' in seconds, a list of time
        parameters which gives the beginning of each frequency change 't_new_start_EIS' in seconds, the final time
        'tf_EIS' in seconds, a list of time parameters which gives the estimated time for reaching equilibrium at each
        frequency 'delta_t_break_EIS' in seconds, and a list of time parameters which gives the estimated time for
        measuring the voltage response at each frequency 'delta_t_measurement_EIS' in seconds.
    current_density : function, optional.
        Current density function.
    """

    # Setting the parameters of the step current density function
    delta_t_ini_step = 120*60 # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_step = 30 # (s). Loading time for the step current density function, from 0 to i_step.
    delta_t_break_step = 15 * 60  # (s). Time at i_step current density for the stabilisation of the internal states.
    i_step = 1.5e4 # (A.m-2). Current density for the step current density function.
    step_current_parameters = {'delta_t_ini_step': delta_t_ini_step, 'delta_t_load_step': delta_t_load_step,
                               'delta_t_break_step': delta_t_break_step,'i_step': i_step}

    # Setting the parameters of the polarization current density function
    delta_t_ini_pola = 120 * 60 # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_pola = 30 # (s). Loading time for one step current of the polarisation current density function.
    delta_t_break_pola = 15 * 60 # (s). Breaking time for one step current, for the stabilisation of the internal states.
    delta_i_pola = 0.05e4 # (A.m-2). Current density step for the polarisation current density function.
    pola_current_parameters = {'delta_t_ini_pola': delta_t_ini_pola, 'delta_t_load_pola': delta_t_load_pola,
                               'delta_t_break_pola': delta_t_break_pola, 'delta_i_pola': delta_i_pola}

    # Setting the parameters of the polarization for calibration current density function
    delta_t_ini_pola_cali = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_pola_cali = 30  # (s). Loading time for one step current of the polarisation current density function.
    delta_t_break_pola_cali = 10 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
    pola_current_for_cali_parameters = {'delta_t_ini_pola_cali': delta_t_ini_pola_cali,
                                        'delta_t_load_pola_cali': delta_t_load_pola_cali,
                                        'delta_t_break_pola_cali': delta_t_break_pola_cali}

    # Setting the parameters of the EIS current density function
    i_EIS, ratio_EIS = 1.0e4, 5/100  # (A/m², ). Parameters for the EIS curve.
    f_EIS = -3, 5, 90, 50 # Frequency parameters for the EIS_current density function.
    t_EIS = EIS_parameters(f_EIS)  # Time parameters for the EIS_current density function.

    # Setting the current density function:
    if type_current == "step": current_density = step_current
    elif type_current == "polarization": current_density = polarization_current
    elif type_current == "polarization_for_cali": current_density = polarization_current_for_calibration
    elif type_current == "EIS": current_density = EIS_current
    elif type_current is None: current_density = None  # No current density function is set.
    else: raise ValueError('You have to specify a type_current which is on the list.')

    return (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters,
            i_EIS, ratio_EIS, f_EIS, t_EIS, current_density)


def operating_inputs_function(pola_current_parameters, type_fuel_cell):
    """This function is used to set the operating inputs of the fuel cell system.

    Parameters
    ----------
    pola_current_parameters : dict
        Parameters for the polarization current density function.
    type_fuel_cell : str
        Type of fuel cell system.

    Returns
    -------
    T_des : float
        Desired fuel cell temperature in Kelvin.
    Pa_des : float
        Desired anode pressure in Pascal.
    Pc_des : float
        Desired cathode pressure in Pascal.
    Sa : float
        Stoichiometric ratio of hydrogen.
    Sc : float
        Stoichiometric ratio of oxygen.
    Phi_a_des : float
        Desired anode relative humidity.
    Phi_c_des : float
        Desired cathode relative humidity.
    pola_current_parameters : dict
        Parameters for the polarization current density. It is a dictionary containing:
        - 'delta_t_ini_pola': the initial time (in seconds) at zero current density for the stabilisation of the
        internal states,
        - 'delta_t_load_pola': the loading time (in seconds) for one step current of the polarisation current
        density function,
        - 'delta_t_break_pola': the breaking time (in seconds) for one step current, for the stabilisation of the
        internal states,
        - 'delta_i_pola': the current density step (in A.m-2) for the polarisation current density function.
        - 'i_max_pola': the maximum current density (in A.m-2) for the polarization curve.
    """

    if type_fuel_cell == "manual_setup": # Setup which are not stored in "stored_operating_inputs".
        T_des = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        y_H2_in = 1 # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        i_max_pola = 3.0e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell is None:
        T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, i_max_pola = None, None, None, None, None, None, None, None, None
    else: # Stored setup in "stored_operating_inputs".
        T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, i_max_pola = stored_operating_inputs(type_fuel_cell)

    pola_current_parameters['i_max_pola'] = i_max_pola  # Update the maximum current density for the polarization curve.
    return T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, pola_current_parameters


def physical_parameters(type_fuel_cell):
    """This function is used to set the physical parameters of the fuel cell system.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell system. It can be "EH-31_1.5", "EH-31_2.0", "EH-31_2.25", "EH-31_2.5", "LF",
        or "manual_setup".

    Returns
    -------
    Hacl : float
        Thickness of the anode catalyst layer in meters.
    Hacl : float
        Thickness of the cathode catalyst layer in meters.
    epsilon_mc : float
        Volume fraction of ionomer in the catalyst layer.
    Hmem : float
        Thickness of the membrane in meters.
    Hgdl : float
        Thickness of the gas diffusion layer in meters.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_c : float
        Compression ratio of the GDL.
    Hmpl : float
        Thickness of the microporous layer in meters.
    epsilon_mpl : float
        Porosity of the microporous layer.
    Hagc : float
        Thickness of the anode gas channel in meters.
    Hcgc : float
        Thickness of the cathode gas channel in meters.
    Wagc : float
        Width of the anode gas channel in meters.
    Wcgc : float
        Width of the cathode gas channel in meters.
    Lgc : float
        Length of the gas channel in meters.
    Aact : float
        Active area of the catalyst layer in meters squared.
    e : float
        Capillary exponent.
    i0_c_ref : float
        Reference exchange current density at the cathode in A.m-2.
    kappa_co : float
        Crossover correction coefficient in mol.m-1.s-1.Pa-1.
    kappa_c : float
        Overpotential correction exponent.
    a_slim : float
        One of the limit liquid saturation coefficients: the slop of slim function.
    b_slim : float
        One of the limit liquid saturation coefficients: the intercept of slim function.
    a_switch : float
        One of the limit liquid saturation coefficients: the slop of s_switch function.
    C_dl : float
        Volumetric double layer capacitance in F.m-3.
    """

    if type_fuel_cell == "manual_setup": # Setup which are not stored in "stored_physical_parameters".
        # Fuel cell physical parameters: 𝜔 (which are not controllable by the system)
        #   Catalyst layer
        Aact = 8.5e-3  # m². It is the active area of the catalyst layer.
        Hacl = 8.089e-6  # m. It is the thickness of the anode catalyst layer.
        Hccl = Hacl  # m. It is the thickness of the cathode catalyst layer.
        epsilon_cl = 0.25  # It is the porosity of the catalyst layer, without units.
        epsilon_mc = 0.3949198274842546  # It is the volume fraction of ionomer in the CL.
        #   Membrane
        Hmem = 2e-5  # m. It is the thickness of the membrane.
        #   Gas diffusion layer
        Hgdl = 2e-4  # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl = 0.7011156494971454  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.27052745219052654  # It is the compression ratio of the GDL.
        #   Microporous layer
        Hmpl = 3e-5  # m. It is the thickness of the microporous layer.
        epsilon_mpl = 0.4  # It is the porosity of the microporous layer.
        #   Gas channel
        Hagc = 5e-4  # m. It is the thickness of the anode gas channel.
        Hcgc = Hagc  # m. It is the thickness of the cathode gas channel.
        Wagc = 4.5e-4  # m. It is the width of the anode gas channel.
        Wcgc = Wagc  # m. It is the width of the cathode gas channel.
        Lgc = 9.67  # m. It is the length of the gas channel.
        #   Interaction parameters between water and PEMFC structure
        e = 5.0  # It is the capillary exponent
        #   Voltage polarization
        i0_c_ref = 2.787917581303015  # A.m-2. It is the reference exchange current density at the cathode.
        kappa_co = 29.793535549174077  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 1.6136446641573106  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.0555312850726664, 0.10514269908118055, 0.6365424991141914  # It is the limit
        #                                                               liquid saturation coefficients.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
    else: # Stored setup in "stored_physical_parameters".
        (Hacl, Hccl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_cl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc,
         Wcgc, Lgc, Aact, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl) = stored_physical_parameters(type_fuel_cell)

    return (Hacl, Hccl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_cl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc,
            Wcgc, Lgc, Aact, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl)


def computing_parameters(step_current_parameters, Hgdl, Hmpl, Hacl, type_fuel_cell):
    """This function is used to set the computing parameters of the fuel cell system.

    Parameters
    ----------
    step_current_parameters : dict
        Parameters for the step current density function.
    Hgdl : float
        Thickness of the gas diffusion layer in meters.
    Hacl : float
        Thickness of the anode catalyst layer in meters.
    type_fuel_cell : str
        Type of fuel cell system. It can be "EH-31_1.5", "EH-31_2.0", "EH-31_2.25", "EH-31_2.5", "LF",
        or "manual_setup".

    Returns
    -------
    n_gdl : int
        Number of model nodes placed inside each GDL.
    n_mpl : int
        Number of model nodes placed inside each MPL.
    t_purge : tuple
        Time parameters for purging the system.
        It is a tuple containing the purge time 'purge_time' in seconds, and the time between two purges
        'delta_purge' in seconds.
        delta_t_dyn_step : float
    step_current_parameters : dict
        Parameters for the step current density. It is a dictionary containing:
        - 'delta_t_ini_step': the initial time (in seconds) at zero current density for the stabilisation of the
        internal states,
        - 'delta_t_load_step': the loading time (in seconds) for the step current density function, from 0 to
        i_step,
        - 'delta_t_break_step': the time (in seconds) at i_step current density for the stabilisation of the
        internal states,
        - 'i_step': the current density (in A.m-2) for the step current density function,
        - 'delta_t_dyn_step': the time (in seconds) for dynamic display of the step current density function.
    """

    n_gdl = max(1, int(Hgdl / Hacl / 4))  # It is the number of model points placed inside each GDL.
    #                                       A good compromise is int(Hgdl/Hacl/4).
    n_mpl = max(1, int(Hmpl / Hacl / 2))  # It is the number of model points placed inside each MPL.
    #                                       A good compromise is int(Hmpl/Hacl/1.5).
    t_purge = 0.6, 15  # (s, s). It is the time parameters for purging the system.
    delta_t_dyn_step = 5*60  # (s). Time for dynamic display of the step current density function.

    if type_fuel_cell == "ZSW-GenStack":
        rtol = 1e-5 # Relative tolerance for the system of ODEs solver.
        atol = 1e-8 # Absolute tolerance for the system of ODEs solver.
    elif type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
        rtol = 1e-5  # Relative tolerance for the system of ODEs solver.
        atol = 1e-8  # Absolute tolerance for the system of ODEs solver.
    else:
        raise ValueError('the type_input given is not valid.')

    step_current_parameters['delta_t_dyn_step'] = delta_t_dyn_step # Update the step current parameters.
    return n_gdl, n_mpl, t_purge, rtol, atol, step_current_parameters

# ____________________________________________Unchanged Physical parameters_____________________________________________
""" These parameters remain unchanged no matter the setting configurations."""

# Physical constants
F = 96485.3321233  # C.mol-1. It is the Faraday constant.
R = 8.31446261815324  # J.mol-1.K-1. It is the universal gas constant.
M_O2 = 3.2e-2  # kg.mol-1. It is the molar mass of O2.
M_H2 = 2e-3  # kg.mol-1. It is the molar mass of H2.
M_N2 = 2.8e-2  # kg.mol-1. It is the molar mass of N2.
M_H2O = M_H2 + M_O2 / 2  # kg.mol-1. It is the molar mass of H2O.
gamma = 1.401  # . It is the heat capacity ratio of dry air at 100°C.
gamma_H2 = 1.404  # . It is the heat capacity ratio of H2 at 100°C.

# External environmental parameters
Text = 298  # K. It is the outside temperature.
Pext = 101325  # Pa. It is the outside pressure.
Phi_ext = 0.4  # It is the outside relative humidity.
y_O2_ext = 0.2095  # . It is the molar fraction of O2 in dry air.

# Model parameters for the cell
rho_mem = 1980  # kg.m-3. It is the density of the dry membrane.
M_eq = 1.1  # kg.mol-1. It is the equivalent molar mass of ionomer.
tau_mpl = 2  # It is the pore structure coefficient in the MPL, without units [Gen Inoue 2016 Journal Power Sources].
tau_cl = 4  # It is the pore structure coefficient in the CL, without units [Gen Inoue 2016 Journal Power Sources].
r_s_gdl = 2.0 # It is the exponent pore blockage in the GDL.
r_s_mpl = 2.5 # It is the exponent pore blockage in the MPL.
r_s_cl = 2.5 # It is the exponent pore blockage in the CL.
Dp_gdl = 33.2e-6  # m. It is the pore diameter of the GDL [ZSW GenStack].
Dp_mpl = 17.4e-6 # m. It is the pore diameter of the MPL [morganUnderstandingGasDiffusion2014].
Dp_cl = 0.15e-6  # m. It is the pore diameter of the CL [Ali Malekian 2019 International Journal of Hydrogen Energy].
theta_c_gdl = 120 * math.pi / 180  # radian. It is the contact angle of GDL for liquid water.
theta_c_mpl = 135 * math.pi / 180  # radian. It is the contact angle of MPL for liquid water.
theta_c_cl = 95 * math.pi / 180  # radian. It is the contact angle of CL for liquid water.
gamma_cond = 5e3  # s-1. It is the overall condensation rate constant for water.
gamma_evap = 1e-4  # Pa-1.s-1. It is the overall evaporation rate constant for water.
epsilon_p = 0.11 #. It is the percolation threshold porosity of the GDL.
alpha_p = 0.785 #. It is a fitted value for the effective matter transfer in the GDL, for through plane direction.
Kshape = 2  # . Mathematical factor governing lambda_eq smoothing.

# Model parameters for the voltage calculation
Re = 1e-06  # ohm.m². It is the estimated electron conduction resistance of the circuit.
C_O2ref = 3.39  # mol.m-3. It is the reference concentration of oxygen.
alpha_c = 0.5  # It is the transfer coefficient of the cathode.
E0 = 1.229  # V. It is the standard-state reversible voltage.
Pref = 1e5  # Pa. It is the reference pressure.
Eact = 73.2e3  # J.mol-1. It is the activation energy.

# Model parameters for the heat transfer calculation
#   Thermal conductivities
k_th_gdl = 1.6 # J.m-1.s-1.K-1. It is the thermal conductivity of the GDLs (non-effective ?) [vetterFreeOpenReference2019].
k_th_mpl = 0.27 # J.m-1.s-1.K-1. It is the thermal conductivity of the MPLs (non-effective ?) [kotakaImpactInterfacialWater2014].
k_th_cl = 0.27 # J.m-1.s-1.K-1. It is the thermal conductivity of the CLs (non-effective ?) [vetterFreeOpenReference2019].
k_th_mem = 0.3 # J.m-1.s-1.K-1. It is the thermal conductivity of the membrane (non-effective ?) [vetterFreeOpenReference2019].
#   Specific heat capacities
Cp_gdl = 568 # J.kg-1.K-1. It is the specific heat capacities of the GDLs [wangQuasi2DTransientModel2018].
Cp_mpl = 568 # J.kg-1.K-1. It is the specific heat capacities of the MPLs [yangEffectsOperatingConditions2019].
Cp_cl = 3300 # J.kg-1.K-1. It is the specific heat capacities the CLs [wangQuasi2DTransientModel2018].
Cp_mem = 833 # J.kg-1.K-1. It is the specific heat capacities of the membrane [wangQuasi2DTransientModel2018].
#   Densities
rho_gdl = 1000 # kg.m-3. It is the density of the GDLs [wangQuasi2DTransientModel2018].
rho_mpl = 1000 # kg.m-3. It is the density of the MPLs [yangEffectsOperatingConditions2019].
rho_cl = 1000 # kg.m-3. It is the density of the CLs [wangQuasi2DTransientModel2018].
#   Electrical conductivities
sigma_e_gdl = 1250 # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
sigma_e_mpl = 5000 # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [yangEffectsOperatingConditions2019].
sigma_e_cl = 350 # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
#   Molar entropy of reactions
delta_s_HOR = 0.104  # J.mol-1.K-1. It is the HOR molar reaction entropy [vetterFreeOpenReference2019].
delta_s_ORR = -163.3  # J.mol-1.K-1. It is the ORR molar reaction entropy [vetterFreeOpenReference2019].

# Model parameters for the balance of plant
#   Physical parameters
n_cell = 1 # . It is the number of cell in the stack.
Vsm = 7.0e-3  # m3. It is the supply manifold volume.
Vem = 2.4e-3  # m-3. It is the exhaust manifold volume.
A_T = 1.18e-3  # m². It is the exhaust manifold throttle area
#   Model parameters
tau_cp = 1  # s. It is the air compressor time constant.
tau_hum = 5  # s. It is the humidifier time constant.
Kp = 5e-8  # m².s-1.Pa-1. It is the proportional constant of the PD controller at the back pressure valve.
Kd = 1e-8  # m².Pa-1. It is the derivative constant of the PD controller at the back pressure valve.
C_D = 5e-2  # . It is the throttle discharge coefficient.
Ksm_in = 1.0e-5  # kg.s-1.Pa-1. It is the supply manifold inlet orifice constant.
Ksm_out = 8.0e-6  # kg.s-1.Pa-1. It is the supply manifold outlet orifice constant.
Kem_in = Ksm_out  # kg.s-1.Pa-1. It is the exhaust manifold inlet orifice constant.
Kem_out = Ksm_in  # kg.s-1.Pa-1. It is the exhaust manifold outlet orifice constant.