# -*- coding: utf-8 -*-

"""This file is used to set the parameters of the fuel cell system.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing functions
from configuration.current_densities import step_current, polarization_current, EIS_current
from modules.settings_modules import stored_operating_inputs, stored_physical_parameters, EIS_parameters

# _______________________________________________________Settings_______________________________________________________
def current_density_parameters(type_current):
    """This function is used to set the parameters of the current density which is imposed to the fuel cell system.

    Parameters
    ----------
    type_current : str
        Type of current density which is imposed to the fuel cell system. It can be "step", "polarization" or "EIS".

    Returns
    -------
    t_step : tuple
        Time parameters for the step_current density function.
        It is a tuple containing the initial time 't0_step' in seconds, final time 'tf_step' in seconds, loading time
        'delta_t_load_step' in seconds, and time for dynamic display 'delta_t_dyn_step' in seconds.
    i_step : tuple
        Current parameters for the step_current density function.
        It is a tuple containing the initial and final current density values 'i_ini_step' and 'i_final_step', both in
        A.m-2.
    delta_pola : tuple
        Parameters for the polarization curve.
        It is a tuple containing the loading time 'delta_t_load_pola' in seconds, the breaking time 'delta_t_break_pola'
        in seconds, the current density step 'delta_i_pola' in A.m-2, and the initial breaking time 'delta_t_ini_pola'
        in seconds.
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
    current_density : function
        Current density function.
    """

    t_step = 0, 1000, 50, 10  # (s, s, s, s). Time parameters for the step_current density function.
    i_step = 0.5e4, 1.5e4  # (A.m-2, A.m-2). Current parameters for the step_current density function.
    delta_pola = 30, 30, 0.1e4, 1 * 60  # (s, s, A.m-2, s). Parameters for the polarization curve.
    i_EIS, ratio_EIS = 1.0e4, 5/100  # (A/m², ). Parameters for the EIS curve.
    f_EIS = -3, 5, 90, 50 # Frequency parameters for the EIS_current density function.
    t_EIS = EIS_parameters(f_EIS)  # Time parameters for the EIS_current density function.
    if type_current == "step": current_density = step_current  # It is the current density function.
    elif type_current == "polarization": current_density = polarization_current  # It is the current density function.
    elif type_current == "EIS": current_density = EIS_current  # It is the current density function.
    else: raise ValueError('You have to specify a type_current which is on the list.')

    return t_step, i_step, delta_pola, i_EIS, ratio_EIS, f_EIS, t_EIS, current_density


def operating_inputs(type_fuel_cell):
    """This function is used to set the operating inputs of the fuel cell system.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell system. It can be "EH-31_1.5", "EH-31_2.0", "EH-31_2.25", "EH-31_2.5", "LF",
        or "manual_setup".

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
    i_max_pola : float
        Maximum current density for the polarization curve.
    """

    if type_fuel_cell == "manual_setup": # Setup which are not stored in "stored_operating_inputs".
        T_des = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        i_max_pola = 3.0e4  # A.m-2. It is the maximum current density for the polarization curve.
    else: # Stored setup in "stored_operating_inputs".
        T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, i_max_pola = stored_operating_inputs(type_fuel_cell)

    return T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, i_max_pola


def physical_parameters(type_fuel_cell):
    """This function is used to set the physical parameters of the fuel cell system.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell system. It can be "EH-31_1.5", "EH-31_2.0", "EH-31_2.25", "EH-31_2.5", "LF",
        or "manual_setup".

    Returns
    -------
    Hcl : float
        Thickness of the anode or cathode catalyst layer in meters.
    epsilon_mc : float
        Volume fraction of ionomer in the catalyst layer.
    tau : float
        Pore structure coefficient in the CL.
    Hmem : float
        Thickness of the membrane in meters.
    Hgdl : float
        Thickness of the gas diffusion layer in meters.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_c : float
        Compression ratio of the GDL.
    Hgc : float
        Thickness of the gas channel in meters.
    Wgc : float
        Width of the gas channel in meters.
    Lgc : float
        Length of the gas channel in meters.
    Aact : float
        Active area of the catalyst layer in meters squared.
    e : float
        Capillary exponent.
    Re : float
        Electron conduction resistance of the circuit in ohm.m².
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
        Hcl = 1e-5  # m. It is the thickness of the anode or cathode catalyst layer.
        epsilon_mc = 0.3949198274842546  # It is the volume fraction of ionomer in the CL.
        tau = 1.015639135686993  # It is the pore structure coefficient in the CL, without units.
        #   Membrane
        Hmem = 2e-5  # m. It is the thickness of the membrane.
        #   Gas diffusion layer
        Hgdl = 2e-4  # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl = 0.7011156494971454  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.27052745219052654  # It is the compression ratio of the GDL.
        #   Gas channel
        Hgc = 5e-4  # m. It is the thickness of the gas channel.
        Wgc = 4.5e-4  # m. It is the width of the gas channel.
        Lgc = 9.67  # m. It is the length of the gas channel.
        #   Interaction parameters between water and PEMFC structure
        e = 5.0  # It is the capillary exponent
        #   Voltage polarization
        Re = 5.694464714060734e-07  # ohm.m². It is the electron conduction resistance of the circuit.
        i0_c_ref = 2.787917581303015  # A.m-2. It is the reference exchange current density at the cathode.
        kappa_co = 29.793535549174077  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 1.6136446641573106  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.0555312850726664, 0.10514269908118055, 0.6365424991141914  # It is the limit
        #                                                               liquid saturation coefficients.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
    else: # Stored setup in "stored_physical_parameters".
        (Hcl, epsilon_mc, tau, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, Re, i0_c_ref, kappa_co,
         kappa_c, a_slim, b_slim, a_switch, C_scl) = stored_physical_parameters(type_fuel_cell)

    return (Hcl, epsilon_mc, tau, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, Re, i0_c_ref, kappa_co,
            kappa_c, a_slim, b_slim, a_switch, C_scl)


def computing_parameters(type_current, Hgdl, Hcl):
    """This function is used to set the computing parameters of the fuel cell system.

    Parameters
    ----------
    type_current : str
        Type of current density which is imposed to the fuel cell system. It can be "step", "polarization" or "EIS".
    Hgdl : float
        Thickness of the gas diffusion layer in meters.
    Hcl : float
        Thickness of the anode or cathode catalyst layer in meters.

    Returns
    -------
    max_step : float
        Maximum time step for the resolution of the system of differential equations.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    t_purge : tuple
        Time parameters for purging the system.
        It is a tuple containing the purge time 'purge_time' in seconds, and the time between two purges
        'delta_purge' in seconds.
    """

    if type_current == "polarization":
        max_step = 0.1  # Once the undetermined parameters are at their optimal values, this max_step can be used.
    elif type_current == "step":
        max_step = 0.1  # it is good enough for having graphs without instabilities.
    else:
        max_step = 0.1
    n_gdl = int(Hgdl / Hcl / 4)  # It is the number of model points placed inside each GDL.
    #                              A good value is int(Hgdl/Hcl/4), which is usually around 5.
    t_purge = 0.6, 15  # (s, s). It is the time parameters for purging the system.

    return max_step, n_gdl, t_purge

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
yO2_ext = 0.2095  # . It is the molar fraction of O2 in dry air.

# Model parameters for the cell
rho_mem = 1980  # kg.m-3. It is the density of the dry membrane.
M_eq = 1.1  # kg.mol-1. It is the equivalent molar mass of ionomer.
epsilon_cl = 0.25  # It is the porosity of the catalyst layer, without units.
theta_c_gdl = 120 * math.pi / 180  # radian. It is the contact angle of GDL for liquid water.
theta_c_cl = 95 * math.pi / 180  # radian. It is the contact angle of CL for liquid water.
gamma_cond = 5e3  # s-1. It is the overall condensation rate constant for water.
gamma_evap = 1e-4  # Pa-1.s-1. It is the overall evaporation rate constant for water.
epsilon_p = 0.11 #. It is the percolation threshold porosity of the GDL.
alpha_p = 0.785 #. It is a fitted value for the effective matter transfer in the GDL, for through plane direction.
Kshape = 2  # . Mathematical factor governing lambda_eq smoothing.

# Model parameters for the voltage calculation
C_O2ref = 3.39  # mol.m-3. It is the reference concentration of oxygen.
alpha_c = 0.5  # It is the transfer coefficient of the cathode.
E0 = 1.229  # V. It is the standard-state reversible voltage.
Pref = 1e5  # Pa. It is the reference pressure.
Eact = 73.2e3  # J.mol-1. It is the activation energy.

# Model parameters for the heat transfer calculation
#   Thermal conductivities
k_th_gdl = 1.6 # J.m-1.s-1.K-1. It is the thermal conductivity of the GDLs [vetterFreeOpenReference2019].
k_th_cl = 0.27 # J.m-1.s-1.K-1. It is the thermal conductivity of the CLs [vetterFreeOpenReference2019].
k_th_mem = 0.3 # J.m-1.s-1.K-1. It is the thermal conductivity of the membrane [vetterFreeOpenReference2019].
#   Specific heat capacities
Cp_gdl = 568 # J.kg-1.K-1. It is the specific heat capacities of the GDLs [wangQuasi2DTransientModel2018].
Cp_cl = 3300 # J.kg-1.K-1. It is the specific heat capacities the CLs [wangQuasi2DTransientModel2018].
Cp_mem = 833 # J.kg-1.K-1. It is the specific heat capacities of the membrane [wangQuasi2DTransientModel2018].
#   Densities
rho_gdl = 1000 # kg.m-3. It is the density of the GDLs [wangQuasi2DTransientModel2018].
rho_cl = 1000 # kg.m-3. It is the density of the CLs [wangQuasi2DTransientModel2018].
#   Electrical conductivities
sigma_e_gdl = 1250 # Ω-1.m-1. It is the electrical conductivity of the GDL [vetterFreeOpenReference2019].
sigma_e_cl = 350 # Ω-1.m-1. It is the electrical conductivity of the GDL [vetterFreeOpenReference2019].
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