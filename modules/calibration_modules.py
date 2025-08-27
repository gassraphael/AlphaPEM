# -*- coding: utf-8 -*-

"""This module contains some of the functions required for the parameter calibration.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import os
from colorama import Fore, Style
import numpy as np

# Importing functions
from configuration.current_densities import polarization_current_for_calibration
from calibration.experimental_values import pola_exp_values_calibration


# _________________________________________________Calibration modules__________________________________________________

def parameter_bounds_for_calibration(type_fuel_cell):
    """This function is used to determine the parameter bounds of the fuel cell model for the calibration when a
    registered type_fuel_cell is considered.

       Parameters
       ----------
       type_fuel_cell : str
               Type of fuel cell configuration.

       Returns
       -------
       varbound : list
            List of the bounds on the parameters to calibrate. Each element is a list containing the minimum and
            maximum values of the parameter, and the type of the parameter ('real' or 'int').
       gene_space : list
            List of dictionaries used to define the bounds of the undetermined parameters for pygad. Each
            dictionary contains the 'low' and 'high' values for the parameter, and optionally a 'step' value for
            integer parameters.
    """

    if type_fuel_cell == "ZSW-GenStack":
        #       Fuel cell physical parameters
        Hacl_min, Hacl_max = 6e-6, 1e-5  # m. It is the thickness of the ACL.
        Hccl_min, Hccl_max = 1e-5, 2e-5  # m. It is the thickness of the CCL.
        Hmem_min, Hmem_max = 1e-5, 2e-5  # m. It is the thickness of the membrane.
        epsilon_gdl_min, epsilon_gdl_max = 0.696, 0.880  # It is the anode/cathode GDL porosity, without units.
        epsilon_mpl_min, epsilon_mpl_max = 0.32, 0.54  # It is the anode/cathode MPL porosity, without units.
        epsilon_cl_min, epsilon_cl_max = 0.40, 0.60  # It is the anode/cathode MPL porosity, without units.
        epsilon_mc_min, epsilon_mc_max = 0.40, 0.60  # It is the volume fraction of ionomer in the CL.
        #       Constants based on the interaction between water and the structure
        e_min, e_max = 3, 5  # It is the capillary exponent, and should be an int number.
        #       Voltage polarization
        i0_c_ref_min, i0_c_ref_max = 1e-3, 80  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co_min, kappa_co_max = 0.01, 40  # A.m-2. It is the crossover correction coefficient.
        kappa_c_min, kappa_c_max = 0.25, 4  # It is the overpotential correction exponent.
        #       The bounds on liquid saturation coefficients are constrained to facilitate calibration.
        a_slim_min, a_slim_max = 0.0, 0.2  # It is one of the limit liquid saturation coefficients.
        b_slim_min, b_slim_max = 0.05, 0.4  # It is one of the limit liquid saturation coefficients.
        a_switch_min, a_switch_max = 0.5, 0.95  # It is one of the limit liquid saturation coefficients.
        #       Undetermined parameter which is not considered yet (require the use of EIS curves to be calibrated)
        C_scl_min, C_sl_max = 2e7, 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
        #       Bounds gathering and type
        varbound = [['epsilon_gdl', epsilon_gdl_min, epsilon_gdl_max, 'real'],
                    ['e', e_min, e_max, 'int'],
                    ['i0_c_ref', i0_c_ref_min, i0_c_ref_max, 'real'],
                    ['kappa_co', kappa_co_min, kappa_co_max, 'real'],
                    ['kappa_c', kappa_c_min, kappa_c_max, 'real']]
        gene_space = []  # List used to define the bounds of the undetermined parameters for pygad.
        for i in range(len(varbound)):
            name, min_val, max_val, type_val = varbound[i]
            if type_val == 'int':
                gene_space.append({'low': min_val, 'high': max_val, 'step': 1})
            else:
                gene_space.append({'low': min_val, 'high': max_val})

    elif type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
        #       Fuel cell physical parameters
        Hacl_min, Hacl_max = 8e-6, 2e-5  # m. It is the thickness of the ACL.
        Hmem_min, Hmem_max = 1.5e-5, 5e-5  # m. It is the thickness of the membrane.
        epsilon_gdl_min, epsilon_gdl_max = 0.50, 0.90  # It is the anode/cathode GDL porosity, without units.
        epsilon_mpl_min, epsilon_mpl_max = 0.30, 0.60  # It is the anode/cathode MPL porosity, without units.
        epsilon_cl_min, epsilon_cl_max = 0.12, 0.50  # It is the anode/cathode MPL porosity, without units.
        epsilon_mc_min, epsilon_mc_max = 0.15, 0.40  # It is the volume fraction of ionomer in the CL.
        epsilon_c_min, epsilon_c_max = 0.15, 0.30  # It is the compression ratio of the GDL.
        #       Constants based on the interaction between water and the structure
        e_min, e_max = 3, 5  # It is the capillary exponent, and should be an int number.
        #       Voltage polarization
        i0_c_ref_min, i0_c_ref_max = 1e-3, 80  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co_min, kappa_co_max = 0.01, 40  # A.m-2. It is the crossover correction coefficient.
        kappa_c_min, kappa_c_max = 0.25, 4  # It is the overpotential correction exponent.
        #       The bounds on liquid saturation coefficients are constrained to facilitate calibration.
        a_slim_min, a_slim_max = 0.0, 0.2  # It is one of the limit liquid saturation coefficients.
        b_slim_min, b_slim_max = 0.05, 0.4  # It is one of the limit liquid saturation coefficients.
        a_switch_min, a_switch_max = 0.5, 0.95  # It is one of the limit liquid saturation coefficients.
        #       Undetermined parameter which is not considered yet (require the use of EIS curves to be calibrated)
        C_scl_min, C_sl_max = 2e7, 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
        #       Bounds gathering and type
        varbound = [['Hacl', Hacl_min, Hacl_max, 'real'],
                    ['Hmem', Hmem_min, Hmem_max, 'real'],
                    ['epsilon_gdl', epsilon_gdl_min, epsilon_gdl_max, 'real'],
                    ['epsilon_mc', epsilon_mc_min, epsilon_mc_max, 'real'],
                    ['e', e_min, e_max, 'int'],
                    ['i0_c_ref', i0_c_ref_min, i0_c_ref_max, 'real'],
                    ['kappa_co', kappa_co_min, kappa_co_max, 'real'],
                    ['kappa_c', kappa_c_min, kappa_c_max, 'real']]
        gene_space = []  # List used to define the bounds of the undetermined parameters for pygad.
        for i in range(len(varbound)):
            name, min_val, max_val, type_val = varbound[i]
            if type_val == 'int':
                gene_space.append({'low': min_val, 'high': max_val, 'step': 1})
            else:
                gene_space.append({'low': min_val, 'high': max_val})
    else:
        raise ValueError("A correct type_fuel_cell should be given.")

    return varbound, gene_space

def parameters_for_calibration(type_fuel_cell):
    """This function is used to determine the parameters of the fuel cell model for the calibration when a registered
    type_fuel_cell is considered.

    Parameters
    ----------
    type_fuel_cell : str
            Type of fuel cell configuration.

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
    Aact : float
        Active area of the cell in m².
    Hmem : float
        Thickness of the membrane in m.
    Hacl : float
        Thickness of the anode catalyst layer in m.
    Hccl : float
        Thickness of the cathode catalyst layer in m.
    Hagc : float
        Thickness of the anode gas channel in m.
    Hcgc : float
        Thickness of the cathode gas channel in m.
    Wagc : float
        Width of the gas anode channel in m.
    Wcgc : float
        Width of the gas cathode channel in m.
    Lgc : float
        Length of the gas channel in m.
    type_auxiliary : str
        Type of auxiliary system.
    type_control : str
        Type of control system.
    type_purge : str
        Type of purge system.
    type_display : str
        Type of display.
    type_plot : str
        Type of plot.
    type_current : str
        Type of current density function.
    current_density : function
        Current density evolution over time. It is a function of time and parameters dictionary.
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
        Current for which a ratio_EIS perturbation is added.
    ratio_EIS : float
        Value of the perturbation on the current density for building the EIS curve.
    t_EIS : tuple
        EIS parameters. It is a tuple containing the initial EIS time after stack equilibrium 't0_EIS', a list of time
        parameters which gives the beginning of each frequency change 't_new_start', the final time 'tf_EIS', a list of
        time parameters which gives the estimated time for reaching equilibrium at each frequency 'delta_t_break_EIS',
        and a list of time parameters which gives the estimated time for measuring the voltage response at each
        frequency 'delta_t_measurement_EIS'.
    f_EIS : tuple
        EIS parameters. It is a tuple containing the power of the initial frequency
        'f_power_min': f_min = 10**f_power_min, the power of the final frequency 'f_power_max', the number of
        frequencies tested 'nb_f' and the number of points calculated per specific period 'nb_points'.
    t_purge : tuple
        Time parameters for purging the system.
        It is the purge time interval 'purge_time' and the time between two purges 'delta_purge'.
    n_gdl : int
        Number of points considered in the GDL.
    i_exp : numpy.ndarray
        Experimental values of the current density.
    U_exp : numpy.ndarray
        Experimental values of the voltage.
    """

    if type_fuel_cell == "ZSW-GenStack":
        # Given values by the author
        #       Operating inputs
        T_des = 68 + 273.15  # K. It is the temperature of the fuel cell.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 39.8, 50  # It is the desired relative humidity.
        Pa_des, Pc_des = 2.2e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        #       Fuel cell physical parameters
        Aact = 2.7972e-2  # m². It is the active area of the catalyst layer.
        Hagc = 2.3e-4  # m. It is the thickness of the anode gas channel.
        Hcgc = 3e-4  # m. It is the thickness of the cathode gas channel.
        Wagc = 4.3e-4  # m. It is the width of the anode gas channel.
        Wcgc = 5.32e-4  # m. It is the width of the cathode gas channel.
        Lgc = 23.31  # m. It is the length of the gas channel.
        #       Fuel cell undetermined physical parameters.
        Hgdl = 1.27e-4  # m. It is the thickness of the gas diffusion layer.
        Hmpl = 7e-5  # m. It is the thickness of the microporous layer.
        epsilon_c = 0.2  # It is the compression ratio of the GDL.

        # Estimated undetermined parameters for the initialisation
        #   Gas diffusion layer
        epsilon_gdl = 0.788  # It is the anode/cathode GDL porosity.
        epsilon_mpl = 0.425  # It is the porosity of the microporous layer.
        #   Catalyst layer
        Hacl = 8e-6  # m. It is the thickness of the anode catalyst layer.
        Hccl = 1.7e-5  # m. It is the thickness of the cathode catalyst layer.
        epsilon_cl = 0.5  # It is the porosity of the microporous layer.
        epsilon_mc = 0.5  # It is the volume fraction of ionomer in the CL.
        #   Membrane
        Hmem = 1.5e-5  # m. It is the thickness of the membrane.
        #   Interaction parameters between water and PEMFC structure
        e = 3.0  # It is the capillary exponent
        #   Voltage polarization
        i0_c_ref = 14.86  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co = 1  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 0.6386  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.05553, 0.10514, 0.63654  # It is the limit liquid saturation coefficients.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.

        # Algorithm parameters for polarization curve generation
        type_auxiliary = "forced-convective_cathode_with_flow-through_anode"
        type_control = "no_control"
        type_purge = "no_purge"
        type_display = "no_display"
        type_plot = "fixed"
        type_current = "polarization_for_cali"
        current_density = polarization_current_for_calibration
        delta_t_ini_step = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_step = 1e-15  # (s). Loading time for the step current density function, from 0 to i_step.
        delta_t_break_step = 0  # (s). Time at i_step current density for the stabilisation of the internal states.
        i_step = 0  # (A.m-2). Current density for the step current density function.
        step_current_parameters = {'delta_t_ini_step': delta_t_ini_step, 'delta_t_load_step': delta_t_load_step,
                                   'delta_t_break_step': delta_t_break_step, 'i_step': i_step}
        delta_t_ini_pola = 30 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        delta_i_pola = 0.05e4  # (A.m-2). Current density step for the polarisation current density function.
        i_max_pola = 1.7e4  # (A.m-2). It is the maximum current density for the polarization curve.
        pola_current_parameters = {'delta_t_ini_pola': delta_t_ini_pola, 'delta_t_load_pola': delta_t_load_pola,
                                   'delta_t_break_pola': delta_t_break_pola, 'delta_i_pola': delta_i_pola,
                                   'i_max_pola': i_max_pola}
        delta_t_ini_pola_cali = 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola_cali = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola_cali = 10 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        pola_current_for_cali_parameters = {'delta_t_ini_pola_cali': delta_t_ini_pola_cali,
                                            'delta_t_load_pola_cali': delta_t_load_pola_cali,
                                            'delta_t_break_pola_cali': delta_t_break_pola_cali}
        i_EIS, ratio_EIS = np.nan, np.nan  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
        f_EIS, t_EIS = np.nan, np.nan  # It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        n_gdl = int(Hgdl / Hacl / 4)  # It is the number of model points placed inside each GDL.
        rtol = 1e-7  # Relative tolerance for the system of ODEs solver.
        atol = 1e-11  # Absolute tolerance for the system of ODEs solver.

    elif type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
        # Given values by the author
        #       Operating inputs
        T_des = 74 + 273.15  # K. It is the temperature of the fuel cell.
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        if type_fuel_cell == "EH-31_1.5":
            Pa_des, Pc_des = 1.5e5, 1.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        elif type_fuel_cell == "EH-31_2.0":
            Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        elif type_fuel_cell == "EH-31_2.25":
            Pa_des, Pc_des = 2.25e5, 2.25e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        else:  # type_fuel_cell == "EH-31_2.5":
            Pa_des, Pc_des = 2.5e5, 2.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        #       Fuel cell physical parameters
        Aact = 8.5e-3  # m². It is the active area of the catalyst layer.
        Wagc = 4.5e-4  # m. It is the width of the anode gas channel.
        Wcgc = Wagc  # m. It is the width of the cathode gas channel.
        Lgc = 9.67  # m. It is the length of the gas channel.

        # Extrapolated physical parameters
        Hgdl = 2e-4  # m. It is the thickness of the gas diffusion layer.
        Hmpl = 3e-5  # m. It is the thickness of the microporous layer.
        epsilon_mpl = 0.4  # It is the porosity of the microporous layer.
        Hagc = 5e-4  # m. It is the thickness of the anode gas channel.
        Hcgc = Hagc  # m. It is the thickness of the cathode gas channel.

        # Estimated undetermined parameters for the initialisation
        #   Gas diffusion layer
        epsilon_gdl = 0.7943  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.2  # It is the compression ratio of the GDL.
        #   Catalyst layer
        Hacl = 8.089e-6  # m. It is the thickness of the anode catalyst layer.
        Hccl = Hacl  # m. It is the thickness of the cathode catalyst layer.
        epsilon_cl = 0.25  # It is the porosity of the catalyst layer, without units.
        epsilon_mc = 0.2111  # It is the volume fraction of ionomer in the CL.
        #   Membrane
        Hmem = 1.5e-5  # m. It is the thickness of the membrane.
        #   Interaction parameters between water and PEMFC structure
        e = 3.0  # It is the capillary exponent
        #   Voltage polarization
        i0_c_ref = 14.86  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co = 1  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 0.6386  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.05553, 0.10514, 0.63654  # It is the limit liquid saturation coefficients.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.

        # Algorithm parameters for polarization curve generation
        type_auxiliary = "forced-convective_cathode_with_flow-through_anode"
        type_control = "no_control"
        type_purge = "no_purge"
        type_display = "no_display"
        type_plot = "fixed"
        type_current = "polarization_for_cali"
        current_density = polarization_current_for_calibration
        delta_t_ini_step = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_step = 1e-15  # (s). Loading time for the step current density function, from 0 to i_step.
        delta_t_break_step = 0  # (s). Time at i_step current density for the stabilisation of the internal states.
        i_step = 0  # (A.m-2). Current density for the step current density function.
        step_current_parameters = {'delta_t_ini_step': delta_t_ini_step, 'delta_t_load_step': delta_t_load_step,
                                   'delta_t_break_step': delta_t_break_step, 'i_step': i_step}
        delta_t_ini_pola = 30 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        delta_i_pola = 0.05e4  # (A.m-2). Current density step for the polarisation current density function.
        i_max_pola = 1.7e4  # (A.m-2). It is the maximum current density for the polarization curve.
        pola_current_parameters = {'delta_t_ini_pola': delta_t_ini_pola, 'delta_t_load_pola': delta_t_load_pola,
                                   'delta_t_break_pola': delta_t_break_pola, 'delta_i_pola': delta_i_pola,
                                   'i_max_pola': i_max_pola}
        delta_t_ini_pola_cali = 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola_cali = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola_cali = 10 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        pola_current_for_cali_parameters = {'delta_t_ini_pola_cali': delta_t_ini_pola_cali,
                                            'delta_t_load_pola_cali': delta_t_load_pola_cali,
                                            'delta_t_break_pola_cali': delta_t_break_pola_cali}
        i_EIS, ratio_EIS = np.nan, np.nan  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
        f_EIS, t_EIS = np.nan, np.nan  # It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        n_gdl = int(Hgdl / Hacl / 4)  # It is the number of model points placed inside each GDL.
        rtol = 1e-7  # Relative tolerance for the system of ODEs solver.
        atol = 1e-11  # Absolute tolerance for the system of ODEs solver.

    else:
        ValueError("A correct type_fuel_cell should be given.")

    # Initialize the operating inputs and parameters dictionaries.
    operating_inputs = {'current_density': current_density, 'T_des': T_des, 'Pa_des': Pa_des, 'Pc_des': Pc_des,
                        'Sa': Sa, 'Sc': Sc, 'Phi_a_des': Phi_a_des, 'Phi_c_des': Phi_c_des}
    current_parameters = {'step_current_parameters': step_current_parameters,
                          'pola_current_parameters': pola_current_parameters,
                          'pola_current_for_cali_parameters': pola_current_for_cali_parameters,
                          'i_EIS': i_EIS, 'ratio_EIS': ratio_EIS, 't_EIS': t_EIS, 'f_EIS': f_EIS}
    accessible_physical_parameters = {'Aact': Aact, 'Hagc': Hagc, 'Hcgc': Hcgc, 'Wagc': Wagc, 'Wcgc': Wcgc, 'Lgc': Lgc}
    undetermined_physical_parameters = {'Hgdl': Hgdl, 'Hmpl': Hmpl, 'Hmem': Hmem, 'Hacl': Hacl, 'Hccl': Hccl,
                                        'epsilon_gdl': epsilon_gdl, 'epsilon_cl': epsilon_cl, 'epsilon_mpl': epsilon_mpl,
                                        'epsilon_mc': epsilon_mc, 'epsilon_c': epsilon_c, 'e': e,
                                        'kappa_co': kappa_co, 'i0_c_ref': i0_c_ref, 'kappa_c': kappa_c,
                                        'a_slim': a_slim, 'b_slim': b_slim, 'a_switch': a_switch,
                                        'C_scl': C_scl}
    computing_parameters = {'n_gdl': n_gdl, 't_purge': t_purge, 'rtol': rtol, 'atol': atol,
                            'type_fuel_cell': type_fuel_cell, 'type_current': type_current,
                            'type_auxiliary': type_auxiliary, 'type_control': type_control, 'type_purge': type_purge,
                            'type_display': type_display, 'type_plot': type_plot}

    # Characteristic points of the experimental polarization curve
    i_exp, U_exp = pola_exp_values_calibration(type_fuel_cell)

    return (operating_inputs, current_parameters, accessible_physical_parameters, undetermined_physical_parameters,
            computing_parameters, i_exp, U_exp)

def update_undetermined_parameters(type_fuel_cell, solution, varbound, undetermined_physical_parameters):
    """
    Update the undetermined physical parameters dictionary with values from the solution.

    Parameters
    ----------
    solution : list
        List of parameter values obtained from the optimization algorithm.
    varbound : list
        List of parameter bounds and names. Each element contains the parameter name at index 0.
    undetermined_physical_parameters : dict
        Dictionary of undetermined physical parameters to be updated.

    Returns
    -------
    dict
        Updated dictionary of undetermined physical parameters.
    """
    for i in range(len(solution)):
        param_name = varbound[i][0]
        if param_name in undetermined_physical_parameters:
            undetermined_physical_parameters[param_name] = solution[i]
        if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
            undetermined_physical_parameters['Hccl'] = undetermined_physical_parameters['Hacl']
    return undetermined_physical_parameters

def calculate_simulation_error(Simulator_1, U_exp_1, i_exp_1, Simulator_2, U_exp_2, i_exp_2):
    """This function is used to calculate the simulation maximal error between the experimental and the simulated
    polarization curves. Two simulations on different operating conditions and on the same stack, and so two set of
    experimental data, are considered as it is the minimum amount of data which is required for the calibration.

    Parameters
    ----------
    Simulator_1 : AlphaPEM object
        PEM simulator which contains the simulation results for the first simulation.
    U_exp_1 : numpy.ndarray
        Experimental values of the voltage for the first simulation.
    i_exp_1 : numpy.ndarray
        Experimental values of the current density for the first simulation.
    Simulator_2 : AlphaPEM object
        PEM simulator which contains the simulation results for the second simulation.
    U_exp_2 : numpy.ndarray
        Experimental values of the voltage for the second simulation.
    i_exp_2 : numpy.ndarray
        Experimental values of the current density for the second simulation.

    Returns
    -------
    sim_error : float
        Maximum error between the experimental and the simulated polarization curves in percentage.
    """

    # Recovery of ifc_1
    t1 = np.array(Simulator_1.variables['t'])
    n1 = len(t1)
    ifc_t_1 = np.zeros(n1)
    for i in range(n1):  # Creation of ifc_t
        ifc_t_1[i] = Simulator_1.operating_inputs['current_density'](t1[i], Simulator_1.parameters)
    # Recovery of ifc_2
    t2 = np.array(Simulator_2.variables['t'])
    n2 = len(t2)
    ifc_t_2 = np.zeros(n2)
    for i in range(n2):  # Creation of ifc_t
        ifc_t_2[i] = Simulator_2.operating_inputs['current_density'](t2[i], Simulator_2.parameters)

    # Polarisation curve point recovery after stack stabilisation for Simulator1
    #   Extraction of the parameters
    #       The initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola_cali_1 = Simulator_1.parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']  # (s).
    #       The loading time for one step current of the polarisation current density function.
    delta_t_load_pola_cali_1 = Simulator_1.parameters['pola_current_for_cali_parameters']['delta_t_load_pola_cali']  # (s).
    #       The breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola_cali_1 = Simulator_1.parameters['pola_current_for_cali_parameters']['delta_t_break_pola_cali']  # (s).
    #   Calculation
    nb_loads1 = len(i_exp_1)  # Number of load which are made
    delta_t_cali_1 = delta_t_load_pola_cali_1 + delta_t_break_pola_cali_1  # s. It is the time of one load.
    ifc_discretized1 = np.zeros(nb_loads1)
    Ucell_discretized1 = np.zeros(nb_loads1)
    for i in range(nb_loads1):
        t_load_1 = delta_t_ini_pola_cali_1 + (i + 1) * delta_t_cali_1 # time for measurement
        idx1 = (np.abs(t1 - t_load_1)).argmin()  # the corresponding index
        ifc_discretized1[i] = ifc_t_1[idx1]  # the last value at the end of each load
        Ucell_discretized1[i] = Simulator_1.variables['Ucell'][idx1]  # the last value at the end of each load
    # Polarisation curve point recovery after stack stabilisation for Simulator2
    #   Extraction of the parameters
    #       The initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola_cali_2 = Simulator_2.parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']  # (s).
    #       The loading time for one step current of the polarisation current density function.
    delta_t_load_pola_cali_2 = Simulator_2.parameters['pola_current_for_cali_parameters']['delta_t_load_pola_cali']  # (s).
    #       The breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola_cali_2 = Simulator_2.parameters['pola_current_for_cali_parameters']['delta_t_break_pola_cali']  # (s).
    #   Calculation
    nb_loads2 = len(i_exp_2)  # Number of load which are made
    delta_t_cali_2 = delta_t_load_pola_cali_2 + delta_t_break_pola_cali_2  # s. It is the time of one load.
    ifc_discretized2 = np.zeros(nb_loads2)
    Ucell_discretized2 = np.zeros(nb_loads2)
    for i in range(nb_loads2):
        t_load_2 = delta_t_ini_pola_cali_2 + (i + 1) * delta_t_cali_2 # time for measurement
        idx2 = (np.abs(t2 - t_load_2)).argmin()  # the corresponding index
        ifc_discretized2[i] = ifc_t_2[idx2]  # the last value at the end of each load
        Ucell_discretized2[i] = Simulator_2.variables['Ucell'][idx2]  # the last value at the end of each load

    # Distance between the simulated and the experimental polarization curves.
    sim_error = (np.max(np.abs(Ucell_discretized1 - U_exp_1) / U_exp_1 * 100)
                     + np.max(np.abs(Ucell_discretized2 - U_exp_2) / U_exp_2 * 100)) / 2  # in %.

    return sim_error


def print_calibration_results(convergence, ga_instance, solution, varbound, sim_error):
    """
    This function prints the calibration results by associating each optimized value with its parameter name.

    Parameters
    ----------
    convergence : dict
        Dictionary containing the convergence information of the genetic algorithm.
    ga_instance : PyGAD object
        Instance of PyGAD used for optimization.
    solution : list
        List of optimized parameter values.
    varbound : list
        List of parameter bounds and names.
    sim_error : float
        Maximum simulation error in percentage.
    """
    print("Convergence:\n", convergence)
    for idx, val in enumerate(solution):
        param_name = varbound[idx][0]
        print(f"Optimized parameter {param_name}: {val}")
    print(Fore.RED + "\nMax simulation error: ", sim_error, "%")
    print(Style.RESET_ALL)
    if ga_instance.best_solution_generation != -1:
        print(f"Best fitness value reached after {ga_instance.best_solution_generation} generations.")


def save_calibration_results(convergence, ga_instance, solution, varbound, sim_error, type_fuel_cell):
    """
    This function saves the calibration results in a text file and a PyGAD file.

    The optimized values are retrieved from the solution list and associated with their names via varbound.

    Parameters
    ----------
    convergence : dict
        Convergence information from the genetic algorithm.
    ga_instance : PyGAD object
        Instance of PyGAD used for optimization.
    solution : list
        List of optimized parameter values.
    varbound : list
        List of parameter bounds and names.
    sim_error : float
        Maximum simulation error in percentage.
    type_fuel_cell : str
        Type of fuel cell configuration.

    Returns
    -------
    None
    """
    root_folder, filename = "results", "parameter_calibration_1.txt"
    subfolder_name = type_fuel_cell[:type_fuel_cell.find('_')] if type_fuel_cell.find('_') != -1 else type_fuel_cell
    counter = 1
    folder_name = os.path.join(root_folder, subfolder_name)
    # Create the folder if necessary
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    # Create the file without erasing the previous ones
    while os.path.isfile(os.path.join(folder_name, filename)):
        counter += 1
        filename = "parameter_calibration_" + str(counter) + ".txt"
    file_path = os.path.join(folder_name, filename)
    # Write information
    with open(file_path, "w") as file:
        file.write("Convergence: " + str(convergence))
        for idx, val in enumerate(solution):
            param_name = varbound[idx][0]
            file.write(f"\nOptimized parameter {param_name}: {val}")
        file.write("\nMax simulation error: " + str(sim_error) + "%")
        file.write("\nAlgorithm works with " + type_fuel_cell + ".")
        if ga_instance.best_solution_generation != -1:
            file.write(f"\nBest fitness value reached after {ga_instance.best_solution_generation} generations.")
    ga_instance.save(filename=os.path.join(folder_name, "parameter_calibration_" + str(counter)))
    if os.path.isfile('parameter_calibration_ongoing.pkl'):
        os.remove('parameter_calibration_ongoing.pkl')