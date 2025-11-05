# -*- coding: utf-8 -*-

"""
- Objectives: Create an open-source software package to simulate the PEM fuel cell for control system applications.
- Authors: Raphaël GASS, Zhongliang LI, Rachid OUTBIB, Samir JEMEI and Daniel HISSEL.
---
This file describes the AlphaPEM class, which is a PEM fuel cell system simulator.
The model is one-dimensional, dynamic, biphasic, and isothermal. It has been published in the following articles:
- Gass et al 2024 J. Electrochem. Soc. https://doi.org/10.1149/1945-7111/ad305a
- Gass et al 2024 SSRN http://dx.doi.org/10.2139/ssrn.4812343
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import os
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Importing constants' value and functions
from model.dif_eq import dydt
from model.flows import calculate_flows
from model.cell_voltage import calculate_cell_voltage
from model.control import control_operating_conditions
from configuration.settings import Pext, Phi_ext, y_O2_ext, C_O2ref, alpha_c, F, R
from modules.dif_eq_modules import event_negative
from modules.transitory_functions import lambda_eq, k_H2, k_O2
from modules.display_modules import (plot_ifc, plot_J, plot_C_v, plot_lambda, plot_s, plot_C_O2, plot_C_H2, plot_C_N2,
                                     plot_T, plot_Ucell, plot_P, plot_Phi_a, plot_Phi_c, plot_Phi_des, plot_v, plot_Re_nb,
                                     plot_polarisation_curve, plot_polarisation_curve_for_cali,
                                     make_Fourier_transformation, plot_EIS_curve_Nyquist, plot_EIS_curve_Bode_amplitude,
                                     plot_EIS_curve_Bode_angle, plot_EIS_curve_tests, plot_power_density_curve,
                                     plot_cell_efficiency, plot_f_drop)
from calibration.experimental_values import pola_exp_values_calibration

# _______________________________________________________AlphaPEM_______________________________________________________

class AlphaPEM:
    def __init__(self, operating_inputs, current_parameters, accessible_physical_parameters,
                 undetermined_physical_parameters, computing_parameters, initial_variable_values=None,
                 time_interval=None):
        """Initialise all parameters defining a fuel cell stack operation: nominal operating conditions,
        applied electrical load, dimensions, and undetermined variables.

        Parameters
        ----------
        operating_inputs : dict
            Dictionary containing the operating inputs for the simulation. It contains:
                - current_density : function
                    Current density evolution over time (operating input). It is a function of time and parameters dictionary.
                - T_des : float
                    Desired fuel cell temperature in Kelvin (operating input).
                - Pa_des : float
                    Desired anode pressure in Pascal (operating input).
                - Pc_des : float
                    Desired cathode pressure in Pascal (operating input).
                - Sa : float
                    Stoichiometric ratio of hydrogen (operating input).
                - Sc : float
                    Stoichiometric ratio of oxygen (operating input).
                - Phi_a_des : float
                    Desired anode relative humidity (operating input).
                - Phi_c_des : float
                    Desired cathode relative humidity (operating input).
        current_parameters : dict
            Dictionary containing the current parameters for the simulation. It contains:
                - step_current_parameters : dict
                    Parameters for the step current density. It is a dictionary containing:
                    - 'delta_t_ini_step': the initial time (in seconds) at zero current density for the stabilisation of the
                    internal states,
                    - 'delta_t_load_step': the loading time (in seconds) for the step current density function, from 0 to
                    i_step,
                    - 'delta_t_break_step': the time (in seconds) at i_step current density for the stabilisation of the
                    internal states,
                    - 'i_step': the current density (in A.m-2) for the step current density function,
                    - 'delta_t_dyn_step': the time (in seconds) for dynamic display of the step current density function.
                - pola_current_parameters : dict
                    Parameters for the polarization current density. It is a dictionary containing:
                    - 'delta_t_ini_pola': the initial time (in seconds) at zero current density for the stabilisation of the
                    internal states,
                    - 'delta_t_load_pola': the loading time (in seconds) for one step current of the polarisation current
                    density function,
                    - 'delta_t_break_pola': the breaking time (in seconds) for one step current, for the stabilisation of the
                    internal states,
                    - 'delta_i_pola': the current density step (in A.m-2) for the polarisation current density function.
                    - 'i_max_pola': the maximum current density (in A.m-2) for the polarization curve.
                - pola_current_for_cali_parameters : dict
                    Parameters for the polarization current density for calibration. It is a dictionary containing:
                    - 'delta_t_ini_pola_cali': the initial time (in seconds) at zero current density for the stabilisation of
                    the internal states,
                    - 'delta_t_load_pola_cali': the loading time (in seconds) for one step current of the polarisation current
                    density function,
                    - 'delta_t_break_pola_cali': the breaking time (in seconds) for one step current, for the stabilisation of
                    the internal states.
                - i_EIS : float
                    Current for which a ratio_EIS perturbation is added (current parameter).
                - ratio_EIS : float
                    Value of the perturbation on the current density for building the EIS curve (current parameter).
                - t_EIS : tuple
                    EIS parameters (current parameters). It is a tuple containing the initial EIS time after stack equilibrium
                    't0_EIS', a list of time parameters which gives the beginning of each frequency change 't_new_start_EIS',
                    the final time 'tf_EIS', a list of time parameters which gives the estimated time for reaching equilibrium
                    at each frequency 'delta_t_break_EIS', and a list of time parameters which gives the estimated time for
                    measuring the voltage response at each frequency 'delta_t_measurement_EIS'.
                f_EIS : tuple
                    EIS parameters (current parameters). It is a tuple containing the power of the initial frequency
                    'f_power_min_EIS': f_min_EIS = 10**f_power_min_EIS, the power of the final frequency 'f_power_max_EIS', the
                    number of frequencies tested 'nb_f_EIS' and the number of points calculated per specific period
                    'nb_points_EIS'.
        accessible_physical_parameters : dict
            Dictionary containing the accessible physical parameters for the simulation. It contains:
                - Aact : float
                    Active area of the cell in m² (accessible physical parameter).
                - n_cell : int
                    Number of cells in the stack (accessible physical parameter).
                - Hagc : float
                    Thickness of the anode gas channel in m (accessible physical parameter).
                Hcgc : float
                    Thickness of the cathode gas channel in m (accessible physical parameter).
                Wagc : float
                    Width of the anode gas channel in m (accessible physical parameter).
                Wcgc : float
                    Width of the cathode gas channel in m (accessible physical parameter).
                Lgc : float
                    Length of the gas channel in m (accessible physical parameter).
        undetermined_physical_parameters : dict
            Dictionary containing the undetermined physical parameters for the simulation. It contains:
                - Hgdl : float
                    Thickness of the gas diffusion layer in m (undetermined physical parameter).
                - Hmem : float
                    Thickness of the membrane in m (undetermined physical parameter).
                - Hacl : float
                    Thickness of the anode catalyst layer in m (undetermined physical parameter).
                - Hccl : float
                    Thickness of the cathode catalyst layer in m (undetermined physical parameter).
                - epsilon_gdl : float
                    Anode/cathode GDL porosity (undetermined physical parameter).
                - epsilon_mc : float
                    Volume fraction of ionomer in the CL (undetermined physical parameter).
                - epsilon_c : float
                    Compression ratio of the GDL (undetermined physical parameter).
                - e : float
                    Capillary exponent (undetermined physical parameter).
                - Re : float
                    Electron conduction resistance in Ω.m2 (undetermined physical parameter).
                - i0_d_c_ref : float
                    Dry reference exchange current density at the cathode in A.m-2 (undetermined physical parameter).
                - kappa_co : float
                    Crossover correction coefficient in mol.m-1.s-1.Pa-1 (undetermined physical parameter).
                - kappa_c : float
                    Overpotential correction exponent (undetermined physical parameter).
                - a_slim : float
                    One of the limit liquid saturation coefficients: the slop of slim function
                    (undetermined physical parameter).
                - b_slim : float
                    One of the limit liquid saturation coefficients: the intercept of slim function
                    (undetermined physical parameter).
                - a_switch : float
                    One of the limit liquid saturation coefficients: the slop of s_switch function
                    (undetermined physical parameter).
                - C_scl : float
                    Volumetric space-charge layer capacitance in F.m-3 (undetermined physical parameter).
        computing_parameters : dict
            Dictionary containing the computing parameters for the simulation. It contains:
                - Htl : float
                    Thickness of the transition layers in meters (computing parameter).
                - n_gdl : int
                    Number of points considered in the GDL (computing parameter).
                - n_mpl : int
                    Number of points considered in the MPL (computing parameter).
                - n_tl : int
                    Number of points considered in the transitory layer (computing parameter).
                - t_purge : tuple
                    Time parameters for purging the system (computing parameter).
                    It is the purge time interval 'purge_time' and the time between two purges 'delta_purge'.
                - rtol : float
                    Relative tolerance for the solver (computing parameter).
                - atol : float
                    Absolute tolerance for the solver (computing parameter).
                - type_fuel_cell : str
                    Type of fuel cell configuration (computing parameter).
                - type_current : str
                    Type of current density function (computing parameter).
                - type_auxiliary : str
                    Type of auxiliary system (computing parameter).
                - type_control : str
                    Type of control system (computing parameter).
                - type_purge : str
                    Type of purge system (computing parameter).
                - type_display : str
                    Type of display (computing parameter).
                - type_plot : str
                    Type of plot (computing parameter).
        initial_variable_values : list, optional
            Initial values of the solver variables. The default is None, which implies that initial values are
            generated considering an equilibrium at the operating inputs without current.
        time_interval : list, optional
            Time intervals for numerical resolution. The default is None, which implies that it is automatically
            generated according to the data given in the current density parameters.
        """

        # Initialize the operating inputs and parameters dictionaries.
        self.operating_inputs = operating_inputs
        self.current_parameters = current_parameters
        self.accessible_physical_parameters = accessible_physical_parameters
        self.undetermined_physical_parameters = undetermined_physical_parameters
        self.computing_parameters = computing_parameters
        self.parameters = {**self.current_parameters, **self.accessible_physical_parameters,
                           **self.undetermined_physical_parameters, **self.computing_parameters}
        if self.operating_inputs['Pa_des'] < Pext or self.operating_inputs['Pc_des'] < Pext:
            raise ValueError('The desired pressure is too low. It cannot be lower than the pressure outside the stack.')

        # Initialize the variables' dictionary.
        self.solver_variable_names = ['C_v_agc', 'C_v_agdl', 'C_v_ampl', 'C_v_acl', 'C_v_ccl', 'C_v_cmpl', 'C_v_cgdl',
                                      'C_v_cgc', 's_agdl', 's_ampl', 's_acl', 's_ccl', 's_cmpl', 's_cgdl', 'lambda_acl',
                                      'lambda_mem', 'lambda_ccl', 'C_H2_agc', 'C_H2_agdl', 'C_H2_ampl', 'C_H2_acl',
                                      'C_O2_ccl', 'C_O2_cmpl', 'C_O2_cgdl', 'C_O2_cgc', 'C_N2_agc', 'C_N2_cgc', 'T_agc',
                                      'T_agdl', 'T_ampl', 'T_acl', 'T_mem', 'T_ccl', 'T_cmpl', 'T_cgdl', 'T_cgc',
                                      'eta_c']
        if self.parameters['type_auxiliary'] == "forced-convective_cathode_with_flow-through_anode" or \
                self.parameters['type_auxiliary'] == "forced-convective_cathode_with_anodic_recirculation":
            self.solver_variable_names.extend(['Pasm', 'Paem', 'Pcsm', 'Pcem', 'Phi_asm', 'Phi_aem', 'Phi_csm',
                                               'Phi_cem', 'Wcp', 'Wa_inj', 'Wc_inj', 'Abp_a', 'Abp_c'])

        self.solver_variable_names_extension()  # Several points are considered in each GC, GDL and MPL. This must be
        #                                        inserted into the solver_variable_names.
        self.all_variable_names = self.solver_variable_names + ['t', 'Ucell', 'v_a_in', 'v_c_in', 'Pa_in', 'Pc_in'] + \
                                                               ['Phi_a_des', 'Phi_c_des']
        self.variables = {key: [] for key in self.all_variable_names}

        # Initialize the control_variables dictionary.
        self.control_variables = {'t_control_Phi': 0,
                                  'Phi_a_des': self.operating_inputs['Phi_a_des'],
                                  'Phi_c_des': self.operating_inputs['Phi_c_des']}

        # Temporary action: simulations with auxiliaries are in reconstruction.
        if self.parameters['type_auxiliary'] != "no_auxiliary":
            self.parameters['type_auxiliary'] == "no_auxiliary"
            print("Auxiliary considerations are temporarily removed, as they require reconstruction.")

        # Create the dynamic evolution.
        #       Create time intervals
        if time_interval is not None:  # Initial time interval may be given to the Simulator.
            self.time_interval = time_interval
        else:  # If not, it is automatically generated.
            self.time_interval = self._create_time_interval()

        #       Create the initial variable values
        if initial_variable_values is not None:  # Initial variable values may be given to the Simulator.
            self.initial_variable_values = initial_variable_values
        else:  # If not, they are generated considering an equilibrium at the operating inputs without current.
            self.initial_variable_values = self._create_initial_variable_values()

        #       Resolution of the system of differential equations.
        event_negative.terminal = True  # Integration is stopped if one of the crucial variables becomes negative.
        self.sol = solve_ivp(dydt, self.time_interval, self.initial_variable_values, method='BDF',
                             rtol = self.parameters['rtol'], atol = self.parameters['atol'], events=event_negative,
                             args=(self.operating_inputs, self.parameters, self.solver_variable_names,
                                   self.control_variables))

        #       Recover the variable values calculated by the solver into the dictionary.
        self._recovery()

        #       Calculate the cell voltage after computing the internal states of the cell.
        self.variables["Ucell"].extend(calculate_cell_voltage(self.variables, self.operating_inputs, self.parameters))

    def solver_variable_names_extension(self):
        """Several points are considered in each GDL and must be inserted into the solver_variable_names.
        """

        new_points_location = ['C_v_agc', 'C_v_agdl', 'C_v_ampl', 'C_v_cmpl', 'C_v_cgdl', 'C_v_cgc',
                               's_agdl', 's_ampl', 's_cmpl', 's_cgdl',
                               'C_H2_agc', 'C_H2_agdl', 'C_H2_ampl', 'C_O2_cmpl', 'C_O2_cgdl', 'C_O2_cgc',
                               'C_N2_agc', 'C_N2_cgc',
                               'T_agc', 'T_agdl', 'T_ampl', 'T_cmpl', 'T_cgdl', 'T_cgc']
        for variable in new_points_location:
            index = self.solver_variable_names.index(variable)
            # Delete the previous points
            self.solver_variable_names.pop(index)
            # Increase the number of points
            if variable.endswith('gc'):
                self.solver_variable_names[index:index] = [f'{variable}_{i}' for i in
                                                           range(1, self.parameters['nb_gc'] + 1)]
            if variable.endswith('gdl'):
                self.solver_variable_names[index:index] = [f'{variable}_{i}' for i in
                                                           range(1, self.parameters['nb_gdl'] + 1)]
            elif variable.endswith('mpl'):
                self.solver_variable_names[index:index] = [f'{variable}_{i}' for i in
                                                           range(1, self.parameters['nb_mpl'] + 1)]

    def _create_time_interval(self):
        """Calculate the time intervals for numerical resolution, according to the current chosen,
        if it is not provided.

        Returns
        -------
        list
            Time interval for numerical resolution. It is used when initial_variable_values == None.
        """

        # Extraction of the parameters
        step_current_parameters  = self.parameters['step_current_parameters']
        pola_current_parameters = self.parameters['pola_current_parameters']
        pola_current_for_cali_parameters = self.parameters['pola_current_for_cali_parameters']
        type_fuel_cell, type_current = self.parameters['type_fuel_cell'], self.parameters['type_current']
        voltage_zone = self.parameters['voltage_zone']

        # Recovery of the good time interval
        if type_current == "step":
            t0_interval = 0 # s.
            tf_interval = step_current_parameters['delta_t_ini_step'] + step_current_parameters['delta_t_load_step'] + \
                          step_current_parameters['delta_t_break_step'] # s.
        elif type_current == "polarization":
            # Extraction of the parameters
            delta_t_ini_pola = pola_current_parameters['delta_t_ini_pola']  # (s).
            delta_t_load_pola = pola_current_parameters['delta_t_load_pola']  # (s).
            delta_t_break_pola = pola_current_parameters['delta_t_break_pola']  # (s).
            delta_i_pola = pola_current_parameters['delta_i_pola']  # (A.m-2).Wccp_des,
            i_max_pola = pola_current_parameters['i_max_pola']  # (A.m-2).
            # Calculation
            t0_interval = 0 # s.
            tf_interval = delta_t_ini_pola + int(i_max_pola / delta_i_pola) * (delta_t_load_pola + delta_t_break_pola)
        elif type_current == "polarization_for_cali":
            # Extraction of the parameters
            delta_t_ini_pola_cali = pola_current_for_cali_parameters['delta_t_ini_pola_cali']  # (s).
            delta_t_load_pola_cali = pola_current_for_cali_parameters['delta_t_load_pola_cali']  # (s).
            delta_t_break_pola_cali = pola_current_for_cali_parameters['delta_t_break_pola_cali']  # (s).
            i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(type_fuel_cell, voltage_zone)  # (A.m-2, V).
            # Calculation
            delta_t_cali = delta_t_load_pola_cali + delta_t_break_pola_cali  # s. It is the time of one load.
            t0_interval = 0
            tf_interval = delta_t_ini_pola_cali + len(i_exp_cali_t) * delta_t_cali # s.
        else:  # EIS time_interval is calculated in the main.py file.
            raise ValueError("Please enter a recognized type_current option for calculating the time interval.")

        # To be reviewed
        self.control_variables['t_control_Phi'] = t0_interval

        return [t0_interval, tf_interval]

    def _create_initial_variable_values(self):
        """Create the initial values of the solver variables if it is not provided.
        It is generated considering an equilibrium inside the fuel cell with H2, O2 and N2, at the external pressure,
        humidity and temperature, without flow or current.

        Returns
        -------
        list
            Initial values of the solver variables. It is used when initial_variable_values == None.
        """

        # Extraction of the operating inputs and parameters
        current_density, T_des = self.operating_inputs['current_density'], self.operating_inputs['T_des']
        Pa_des, Pc_des = self.operating_inputs['Pa_des'], self.operating_inputs['Pc_des']
        Phi_a_des, Phi_c_des = self.operating_inputs['Phi_a_des'], self.operating_inputs['Phi_c_des']
        y_H2_in = self.operating_inputs['y_H2_in']
        Hmem, kappa_co, kappa_c = self.parameters['Hmem'], self.parameters['kappa_co'], self.parameters['kappa_c']
        i0_d_c_ref, i0_h_c_ref = self.parameters['i0_d_c_ref'], self.parameters['i0_h_c_ref']
        a_slim, b_slim, a_switch = self.parameters['a_slim'], self.parameters['b_slim'], self.parameters['a_switch']
        nb_gc, nb_gdl, nb_mpl = self.parameters['nb_gc'], self.parameters['nb_gdl'], self.parameters['nb_mpl']
        type_auxiliary = self.parameters['type_auxiliary']

        # Initial fuel cell states
        #   Intermediate values
        T_ini = T_des  # K. It is the initial temperature in the fuel cell. For now, the temperature dynamic evolution
        #                   is not considered inside the gas channel and bipolar plates.
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
                type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            Pa_ini, Pc_ini = Pext, Pext  # Pa. It is the initial pressure in the fuel cell.
            Phi_a_ini, Phi_c_ini = Phi_ext, Phi_ext  # It is the initial relative humidity in the fuel cell.
        else: # type_auxiliary == "no_auxiliaries". For the no_auxiliaries case, the desired values are directly taken as initial values.
            Pa_ini, Pc_ini = Pa_des, Pc_des # Indeed, the dynamics of the auxiliaries are not considered.
            Phi_a_ini, Phi_c_ini = Phi_a_des, Phi_c_des # Choosing other values would create an imbalance in the initial state.
        Psat_ini = 101325 * 10 ** (-2.1794 + 0.02953 * (T_ini - 273.15) - 9.1837e-5 * (T_ini - 273.15) ** 2 +
                                   1.4454e-7 * (T_ini - 273.15) ** 3)
        slim = a_slim * (Pc_ini / 1e5) + b_slim
        s_switch = a_switch * slim
        #   Initial fuel cell states
        C_v_a_ini = Phi_a_ini * Psat_ini / (R * T_ini)  # mol.m-3. It is the initial vapor concentration.
        C_v_c_ini = Phi_c_ini * Psat_ini / (R * T_ini)  # mol.m-3. It is the initial vapor concentration.
        C_O2_ini = y_O2_ext * (Pc_ini - Phi_c_ini * Psat_ini) / (R * T_ini)  # mol.m-3. It is the initial O2
        #                                                                        concentration in the fuel cell.
        C_N2_cgc_ini = (1 - y_O2_ext) * (Pc_ini - Phi_c_ini * Psat_ini) / (R * T_ini)  # mol.m-3. It is the initial N2
        #                                                                          concentration in the anode fuel cell.
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            C_H2_ini = y_H2_in * (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)  # mol.m-3. It is the initial H2
        #                                                                        concentration in the fuel cell.
            C_N2_agc_ini = (1 - y_H2_in) * (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)  # mol.m-3. It is the initial N2
            #                                                                       concentration in the anode fuel cell.
        else:
            C_H2_ini = (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)  # mol.m-3. It is the initial H2
            #                                                                        concentration in the fuel cell.
            C_N2_agc_ini = 0                                            # mol.m-3. It is the initial N2
            #                                                                      concentration in the anode fuel cell.
        s_ini = 0  # It is the initial liquid water saturation in the fuel cell.
        lambda_mem_ini = lambda_eq(C_v_c_ini, s_ini, T_ini)  # It is the initial water content in the fuel cell.
        i_fc_ini = current_density(self.time_interval[0], self.parameters)
        i_n_ini = 2 * F * R * T_ini / Hmem * C_H2_ini * k_H2(lambda_mem_ini, T_ini, kappa_co) + \
                  4 * F * R * T_ini / Hmem * C_O2_ini * k_O2(lambda_mem_ini, T_ini, kappa_co)
        f_drop_ini = 0.5 * (1.0 - math.tanh((4 * s_ini - 2 * slim - 2 * s_switch) / (slim - s_switch)))
        eta_c_ini = R * T_ini / (alpha_c * F) * \
                    math.log((i_fc_ini + i_n_ini) / (i0_d_c_ref ** (1 - f_drop_ini) * i0_h_c_ref ** f_drop_ini) *
                             (C_O2ref / C_O2_ini) ** kappa_c)  # It is the initial
        #                                                                       cathode overpotential in the fuel cell.

        # Initial auxiliary system state
        Wcp_ini = 0  # kg.s-1. It is the flow rate of the air compressor.
        Wa_inj_ini = 0  # kg.s-1. It is the flow rate of the air compressor at the anode side.
        Wc_inj_ini = 0  # kg.s-1. It is the flow rate of the air compressor at the cathode side.
        Abp_a_ini = 0  # It is the throttle area of the back pressure valve at the anode.
        Abp_c_ini = 0  # It is the throttle area of the back pressure valve at the cathode.

        # Main variable initialization
        C_v_agc, C_v_agdl, C_v_ampl, C_v_acl = [C_v_a_ini] * 4
        C_v_ccl, C_v_cmpl, C_v_cgdl, C_v_cgc = [C_v_c_ini] * 4
        s_agdl, s_ampl, s_acl, s_ccl, s_cmpl, s_cgdl = [s_ini] * 6
        s_boundary = 0  # Dirichlet boundary condition
        lambda_acl, lambda_mem, lambda_ccl = [lambda_mem_ini] * 3
        C_H2_agc, C_H2_agdl, C_H2_ampl, C_H2_acl = [C_H2_ini] * 4
        C_O2_ccl, C_O2_cmpl, C_O2_cgdl, C_O2_cgc = [C_O2_ini] * 4
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            C_N2_agc, C_N2_cgc = C_N2_agc_ini, C_N2_cgc_ini # Test bench: simulated H2 recirculation which leads to N2 in the anode.
        else:
            C_N2_agc, C_N2_cgc = 0, C_N2_cgc_ini
        T_agc, T_agdl, T_ampl, T_acl, T_mem, T_ccl, T_cmpl, T_cgdl, T_cgc = [T_ini] * 9
        eta_c = eta_c_ini
        Pasm, Paem = [Pa_ini] * 2
        Pcsm, Pcem = [Pc_ini] * 2
        Phi_asm, Phi_aem = [Phi_a_ini] * 2
        Phi_csm, Phi_cem = [Phi_c_ini] * 2
        Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c = Wcp_ini, Wa_inj_ini, Wc_inj_ini, Abp_a_ini, Abp_c_ini

        # Gathering of the variables initial value into one list
        initial_variable_values = ([C_v_agc] * nb_gc + [C_v_agdl] * nb_gdl + [C_v_ampl] * nb_mpl + [C_v_acl, C_v_ccl] +
                                   [C_v_cmpl] * nb_mpl + [C_v_cgdl] * nb_gdl + [C_v_cgc] * nb_gc +
                                   [s_boundary] + [s_agdl] * (nb_gdl - 1) + [s_ampl] * nb_mpl + [s_acl, s_ccl] +
                                   [s_cmpl] * nb_mpl + [s_cgdl] * (nb_gdl - 1) + [s_boundary] +
                                   [lambda_acl, lambda_mem, lambda_ccl] +
                                   [C_H2_agc] * nb_gc + [C_H2_agdl] * nb_gdl + [C_H2_ampl] * nb_mpl +
                                   [C_H2_acl, C_O2_ccl] + [C_O2_cmpl] * nb_mpl + [C_O2_cgdl] * nb_gdl +
                                   [C_O2_cgc] * nb_gc +
                                   [C_N2_agc] * nb_gc + [C_N2_cgc] * nb_gc +
                                   [T_agc] * nb_gc + [T_agdl] * nb_gdl + [T_ampl] * nb_mpl + [T_acl, T_mem, T_ccl] +
                                   [T_cmpl] * nb_mpl + [T_cgdl] * nb_gdl + [T_cgc] * nb_gc +
                                   [eta_c])
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode" or \
                type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
                initial_variable_values.extend([Pasm, Paem, Pcsm, Pcem] + [Phi_asm, Phi_aem] + [Phi_csm, Phi_cem] +
                                               [Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c])

        return initial_variable_values

    def _recovery(self):
        """Recover the values which have been calculated by the solver and add them into the variables' dictionary.
        However, the numerical resolution method does not, by design, recover all the internal states of the stack,
        even though they are calculated during this process. They therefore have to be recovered manually.
        """

        # Recovery of the time span
        self.variables['t'].extend(list(self.sol.t))

        # Recovery of the main variables dynamic evolution
        for index, key in enumerate(self.solver_variable_names):
            self.variables[key].extend(list(self.sol.y[index]))

        # Recovery of more variables
        if self.parameters['type_display'] != "no_display":
            #   The control variables should be reinitialized. To be reviewed.
            if self.parameters['type_current'] == "step":
                self.control_variables['t_control_Phi'] = 0
            else:
                self.control_variables['t_control_Phi'] = 0
            self.control_variables['Phi_a_des'] = self.operating_inputs['Phi_a_des']
            self.control_variables['Phi_c_des'] = self.operating_inputs['Phi_c_des']

            for j in range(len(self.sol.t)):  # For each time...
                # ... recovery of i_fc.
                i_fc = self.operating_inputs["current_density"](self.variables['t'][j], self.parameters)
                # ... recovery of S_abs_acl, S_abs_ccl, Jmem_acl, Jmem_ccl.
                last_solver_variables = {key: self.variables[key][j] for key in self.solver_variable_names}
                flows_recovery = calculate_flows(self.variables['t'][j], last_solver_variables, self.control_variables,
                                                 i_fc, self.operating_inputs, self.parameters)
                for key in ['v_a_in', 'v_c_in', 'Pa_in', 'Pc_in']:
                    self.variables[key].append(flows_recovery[key])
                # ... recovery of Phi_a_des and Phi_c_des.
                if self.parameters["type_control"] == "Phi_des":
                    sv = {'lambda_mem': self.variables['lambda_mem'][j], 's_ccl': self.variables['s_ccl'][j]}
                    control_operating_conditions(self.variables['t'][j], sv, self.operating_inputs,
                                                 self.parameters, self.control_variables)
                    for key in ['Phi_a_des', 'Phi_c_des']: self.variables[key].append(self.control_variables[key])

    def Display(self, ax1=None, ax2=None, ax3=None):
        """Display the plots of the program.

        Parameters
        ----------
        ax1 : matplotlib.axes.Axes, optional
            Axes for the first set of plots. The default is None.
        ax2 : matplotlib.axes.Axes, optional
            Axes for the second set of plots. The default is None.
        ax3 : matplotlib.axes.Axes, optional
            Axes for the third set of plots. The default is None.
        """

        # Extraction of the operating inputs and parameters
        nb_gdl, type_fuel_cell = self.parameters['nb_gdl'], self.parameters['type_fuel_cell']
        type_current, type_display = self.parameters['type_current'], self.parameters['type_display']

        # Parameters' preparation
        n = len(self.variables['t'])
        subfolder_name = type_fuel_cell[:type_fuel_cell.rfind('_')] if type_fuel_cell.rfind('_') != -1 \
            else type_fuel_cell

        # Display
        if type_current == "step":
            if type_display == "multiple":

                figs, axes = zip(*[plt.subplots(figsize=(8, 8)) for _ in range(13)])

                plot_ifc(self.variables, self.operating_inputs, self.parameters, axes[0])
                plot_J(self.variables, self.parameters, axes[1])
                plot_C_v(self.variables, self.parameters, axes[2])
                plot_lambda(self.variables, self.operating_inputs, self.parameters, axes[3])
                plot_s(self.variables, self.operating_inputs, self.parameters, axes[4])
                plot_C_O2(self.variables, self.parameters, axes[5])
                plot_C_H2(self.variables, self.parameters, axes[6])
                plot_C_N2(self.variables, self.parameters, axes[7])
                plot_T(self.variables, self.operating_inputs, self.parameters, axes[8])
                plot_Ucell(self.variables, self.parameters, axes[9])
                plot_P(self.variables, self.operating_inputs, self.parameters, axes[10])
                plot_Phi_a(self.variables, self.operating_inputs, self.parameters, axes[11])
                plot_Phi_c(self.variables, self.operating_inputs, self.parameters, axes[12])

                # Considering the number of plots, the saving instructions are made here and not in the main.py file.
                self.Saving_instructions("results", subfolder_name, "step_current_ifc_1.pdf", figs[0])
                self.Saving_instructions("results", subfolder_name, "step_current_J_1.pdf", figs[1])
                self.Saving_instructions("results", subfolder_name, "step_current_Cv_1.pdf", figs[2])
                self.Saving_instructions("results", subfolder_name, "step_current_lambda_1.pdf", figs[3])
                self.Saving_instructions("results", subfolder_name, "step_current_s_1.pdf", figs[4])
                self.Saving_instructions("results", subfolder_name, "step_current_C_O2_1.pdf", figs[5])
                self.Saving_instructions("results", subfolder_name, "step_current_C_H2_1.pdf", figs[6])
                self.Saving_instructions("results", subfolder_name, "step_current_C_N2_1.pdf", figs[7])
                self.Saving_instructions("results", subfolder_name, "step_current_T_1.pdf", figs[8])
                self.Saving_instructions("results", subfolder_name, "step_current_Ucell_1.pdf", figs[9])
                self.Saving_instructions("results", subfolder_name, "step_current_P_1.pdf", figs[10])
                self.Saving_instructions("results", subfolder_name, "step_current_Phi_a_1.pdf", figs[11])
                self.Saving_instructions("results", subfolder_name, "step_current_Phi_c_1.pdf", figs[12])

                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "synthetic":

                plot_ifc(self.variables, self.operating_inputs, self.parameters, ax1[0, 0])
                plot_Ucell(self.variables, self.parameters, ax1[0, 1])
                plot_T(self.variables, self.operating_inputs, self.parameters, ax1[0, 2])
                plot_C_v(self.variables, self.parameters, ax1[1, 0])
                plot_s(self.variables, self.operating_inputs, self.parameters, ax1[1, 1])
                plot_lambda(self.variables, self.operating_inputs, self.parameters, ax1[1, 2])
                plot_C_H2(self.variables, self.parameters, ax1[2, 0])
                plot_C_O2(self.variables, self.parameters, ax1[2, 1])
                plot_P(self.variables, self.operating_inputs, self.parameters, ax1[2, 2])

                plot_v(self.variables, self.parameters, ax2[0])
                plot_Re_nb(self.variables, self.parameters, ax2[1])

                plt.pause(1)  # A break is necessary to plot the new points in dynamic mode

        elif type_current == "polarization":
            if type_display == "multiple":

                plot_polarisation_curve(self.variables, self.operating_inputs, self.parameters, ax1[0])
                plot_power_density_curve(self.variables, self.operating_inputs, self.parameters, n, ax1[1])
                plot_cell_efficiency(self.variables, self.operating_inputs, self.parameters, n, ax1[2])

                plot_Phi_des(self.variables, self.operating_inputs, self.parameters, ax2[0])
                plot_lambda(self.variables, self.operating_inputs, self.parameters, ax2[1])
                plot_s(self.variables, self.operating_inputs, self.parameters, ax2[2])
                plot_T(self.variables, self.operating_inputs, self.parameters, ax2[3])

                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "synthetic":

                plot_polarisation_curve(self.variables, self.operating_inputs, self.parameters, ax1[0])
                plot_s(self.variables, self.operating_inputs, self.parameters, ax1[1])
                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "no_display":

                plot_polarisation_curve(self.variables, self.operating_inputs, self.parameters, ax1, show=False)

        elif type_current == "polarization_for_cali":
            if type_display == "multiple":

                plot_polarisation_curve_for_cali(self.variables, self.operating_inputs, self.parameters, ax1[0])
                plot_lambda(self.variables, self.operating_inputs, self.parameters, ax1[1])
                plot_s(self.variables, self.operating_inputs, self.parameters, ax1[2])

                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "synthetic":

                plot_polarisation_curve_for_cali(self.variables, self.operating_inputs, self.parameters, ax1)
                plt.pause(0.01)  # A break is necessary to plot the new points in dynamic mode

        elif type_current == "EIS":
            if type_display == "multiple":

                Fourier_results = make_Fourier_transformation(self.variables, self.operating_inputs, self.parameters)
                plot_EIS_curve_Nyquist(self.parameters, Fourier_results, ax1)
                plot_EIS_curve_Bode_amplitude(self.parameters, Fourier_results, ax2)
                plot_EIS_curve_Bode_angle(self.parameters, Fourier_results, ax3)

                # # Tests to verify the accuracy of EIS simulation.
                # plot_EIS_curve_tests(self.variables, self.operating_inputs, self.parameters, Fourier_results)

                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "synthetic":

                Fourier_results = make_Fourier_transformation(self.variables, self.operating_inputs, self.parameters)
                plot_EIS_curve_Nyquist(self.parameters, Fourier_results, ax1[0])
                plot_EIS_curve_Bode_amplitude(self.parameters, Fourier_results, ax1[1])
                plot_EIS_curve_Bode_angle(self.parameters, Fourier_results, ax1[2])

                # # Tests to verify the accuracy of EIS simulation.
                # plot_EIS_curve_tests(self.variables, self.operating_inputs, self.parameters, Fourier_results)

                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode

    def Save_plot(self, fig1=None, fig2=None, fig3=None):
        """Saves the plots. The names of the files are automatically generated according to the type_current and the
        type_display.

        Parameters
        ----------
        fig1 : matplotlib.figure.Figure, optional
            Figure for the first plot. The default is None.
        fig2 : matplotlib.figure.Figure, optional
            Figure for the second plot. The default is None.
        fig3 : matplotlib.figure.Figure, optional
            Figure for the third plot. The default is None.
        """

        # Extraction of the operating inputs and parameters
        type_fuel_cell, type_current = self.parameters['type_fuel_cell'], self.parameters['type_current']
        type_display = self.parameters['type_display']

        # Folder name
        subfolder_name = type_fuel_cell.split('_')[0]

        # For the step current
        if type_current == "step":
            if type_display == "multiple":
                pass  # saving instruction is directly implemented within AlphaPEM.Display for this situation.
            if type_display == "synthetic":
                self.Saving_instructions("results", subfolder_name, "step_current_syn_1.pdf", fig1)

        # For the polarization curve
        elif type_current == "polarization":
            if type_display == "multiple":
                self.Saving_instructions("results", subfolder_name, "global_indicators_1.pdf", fig1)
                self.Saving_instructions("results", subfolder_name, "pola_curve_syn_1.pdf", fig2)
            elif type_display == "synthetic":
                self.Saving_instructions("results", subfolder_name, "pola_curve_1.pdf", fig1)

        # For the EIS curve
        elif type_current == "EIS":
            if type_display == "multiple":
                self.Saving_instructions("results", subfolder_name, "Nyquist_plot_1.pdf", fig1)
                self.Saving_instructions("results", subfolder_name, "Bode_amplitude_curve_1.pdf", fig2)
                self.Saving_instructions("results", subfolder_name, "Bode_angle_curve_1.pdf", fig3)
            elif type_display == "synthetic":
                self.Saving_instructions("results", subfolder_name, "Nyquist_plot_syn_1.pdf", fig1)

        # For the polarization curve
        elif type_current == "polarization":
            if type_display == "multiple":
                self.Saving_instructions("results", subfolder_name, "impact_cali_on_internal_state_1.pdf", fig1)
            elif type_display == "synthetic":
                self.Saving_instructions("results", subfolder_name, "pola_curve_cali_1.pdf", fig1)

    def Saving_instructions(self, root_folder, subfolder_name, filename, fig):
        """Gives the saving instructions for the figures.

        Parameters
        ----------
        root_folder : str
            The root folder for the saving.
        subfolder_name : str
            The subfolder name for the saving.
        filename : str
            The filename for the saving.
        fig : matplotlib.figure.Figure
            The figure to be saved.
        """

        # Create the folder if necessary
        folder_name = os.path.join(root_folder, subfolder_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # Create the filename without erasing the previous ones
        counter = 1
        while os.path.isfile(os.path.join(folder_name, filename)):
            counter += 1
            if filename[-6] == "_":  # for the numbers between 1 and 9
                filename = filename[:-5] + str(counter) + ".pdf"
            elif filename[-7] == "_":  # for the numbers between 10 and 99.
                filename = filename[:-6] + str(counter) + ".pdf"
            else:  # for the numbers between 100 and 999. The bigger numbers are not considered.
                filename = filename[:-7] + str(counter) + ".pdf"

        # Save the figure
        file_path = os.path.join(folder_name, filename)
        fig.savefig(file_path, dpi=900, transparent=False, bbox_inches='tight')