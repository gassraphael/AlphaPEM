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
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.integrate import solve_ivp

# Importing constants' value and functions
from model.dif_eq import dydt
from model.flows import calculate_flows
from model.cell_voltage import calculate_cell_voltage
from model.control import control_operating_conditions
from configuration.settings import Pext, Kshape, yO2_ext, C_O2ref, alpha_c, F, R
from modules.dif_eq_modules import event_negative
from modules.transitory_functions import lambda_eq, C_v_sat, k_H2, k_O2
from modules.display_modules import plot_ifc, plot_J, plot_C_v, plot_lambda, plot_s, plot_C_O2, plot_C_H2, plot_C_N2, \
    plot_Ucell, plot_P, plot_Phi, plot_Phi_des, plot_polarisation_curve, make_Fourier_transformation, \
    plot_EIS_curve_Nyquist, plot_EIS_curve_Bode_amplitude, plot_EIS_curve_Bode_angle, plot_EIS_curve_tests, \
    plot_power_density_curve, plot_cell_efficiency
from modules.main_modules import saving_instructions

# PyCharm requirement for dynamic plot display
mpl.use("Qt5Agg")


# _______________________________________________________AlphaPEM_______________________________________________________

class AlphaPEM:

    def __init__(self, current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step, i_max_pola,
                 delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau,
                 epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step,
                 n_gdl, t_purge, type_fuel_cell, type_current, type_auxiliary, type_control, type_purge, type_display,
                 type_plot, initial_variable_values=None, time_interval=None):
        """Initialise all parameters defining a fuel cell stack operation: nominal operating conditions,
        applied electrical load, dimensions, and undetermined variables.

        Parameters
        ----------
        current_density : function
            Current density evolution over time (operating input). It is a function of time and parameters dictionary.
        Tfc : float
            Desired fuel cell temperature in Kelvin (operating input).
        Pa_des : float
            Desired anode pressure in Pascal (operating input).
        Pc_des : float
            Desired cathode pressure in Pascal (operating input).
        Sa : float
            Stoichiometric ratio of hydrogen (operating input).
        Sc : float
            Stoichiometric ratio of oxygen (operating input).
        Phi_a_des : float
            Desired anode relative humidity (operating input).
        Phi_c_des : float
            Desired cathode relative humidity (operating input).
        t_step : tuple
            Time parameters for the step_current density function (current parameters).
            It is a tuple containing the initial time 't0_step', final time 'tf_step', loading time 'delta_t_load_step'
            and dynamic time for display 'delta_t_dyn_step'.
        i_step : tuple
            Current parameters for the step_current density function (current parameters).
            It is a tuple containing the initial and final current density value 'i_ini_step' and 'i_final_step'.
        i_max_pola : float
            Maximum current density for the polarization curve (current parameter).
        delta_pola : tuple
            Parameters for the polarization curve (current parameters). It is a tuple containing the loading time
            'delta_t_load_pola', the breaking time 'delta_t_break_pola', the current density step 'delta_i_pola', and
            the initial breaking time 'delta_t_ini_pola'.
        i_EIS : float
            Current for which a ratio_EIS perturbation is added (current parameter).
        ratio_EIS : float
            Value of the perturbation on the current density for building the EIS curve (current parameter).
        t_EIS : tuple
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
        Aact : float
            Active area of the cell in m² (accessible physical parameter).
        Hgdl : float
            Thickness of the gas diffusion layer in m (accessible physical parameter).
        Hmem : float
            Thickness of the membrane in m (accessible physical parameter).
        Hcl : float
            Thickness of the catalyst layer in m (accessible physical parameter).
        Hgc : float
            Thickness of the gas channel in m (accessible physical parameter).
        Wgc : float
            Width of the gas channel in m (accessible physical parameter).
        Lgc : float
            Length of the gas channel in m (accessible physical parameter).
        epsilon_gdl : float
            Anode/cathode GDL porosity (undetermined physical parameter).
        tau : float
            Pore structure coefficient (undetermined physical parameter).
        epsilon_mc : float
            Volume fraction of ionomer in the CL (undetermined physical parameter).
        epsilon_c : float
            Compression ratio of the GDL (undetermined physical parameter).
        e : float
            Capillary exponent (undetermined physical parameter).
        Re : float
            Electron conduction resistance of the circuit in ohm.m² (undetermined physical parameter).
        i0_c_ref : float
            Reference exchange current density at the cathode in A.m-2 (undetermined physical parameter).
        kappa_co : float
            Crossover correction coefficient in mol.m-1.s-1.Pa-1 (undetermined physical parameter).
        kappa_c : float
            Overpotential correction exponent (undetermined physical parameter).
        a_slim : float
            One of the limit liquid saturation coefficients: the slop of slim function
            (undetermined physical parameter).
        b_slim : float
            One of the limit liquid saturation coefficients: the intercept of slim function
            (undetermined physical parameter).
        a_switch : float
            One of the limit liquid saturation coefficients: the slop of s_switch function
            (undetermined physical parameter).
        C_dl : float
            Volumetric double layer capacitance in F.m-3 (undetermined physical parameter).
        max_step : float
            Maximum time step for the solver (computing parameter).
        n_gdl : int
            Number of points considered in the GDL (computing parameter).
        t_purge : tuple
            Time parameters for purging the system (computing parameter).
            It is the purge time interval 'purge_time' and the time between two purges 'delta_purge'.
        type_fuel_cell : str
            Type of fuel cell configuration (computing parameter).
        type_current : str
            Type of current density function (computing parameter).
        type_auxiliary : str
            Type of auxiliary system (computing parameter).
        type_control : str
            Type of control system (computing parameter).
        type_purge : str
            Type of purge system (computing parameter).
        type_display : str
            Type of display (computing parameter).
        type_plot : str
            Type of plot (computing parameter).
        initial_variable_values : list, optional
            Initial values of the solver variables. The default is None, which implies that initial values are
            generated considering an equilibrium at the operating inputs without current.
        time_interval : list, optional
            Time intervals for numerical resolution. The default is None, which implies that it is automatically
            generated according to the data given in the current density parameters.
        """

        # Initialize the operating inputs and parameters dictionaries.
        self.operating_inputs = {'current_density': current_density, 'Tfc': Tfc, 'Pa_des': Pa_des, 'Pc_des': Pc_des,
                                 'Sa': Sa, 'Sc': Sc, 'Phi_a_des': Phi_a_des, 'Phi_c_des': Phi_c_des}
        self.current_parameters = {'t_step': t_step, 'i_step': i_step, 'delta_pola': delta_pola,
                                   'i_max_pola': i_max_pola, 'i_EIS': i_EIS, 'ratio_EIS': ratio_EIS, 't_EIS': t_EIS,
                                   'f_EIS': f_EIS}
        self.accessible_physical_parameters = {'Aact': Aact, 'Hgdl': Hgdl, 'Hmem': Hmem, 'Hcl': Hcl, 'Hgc': Hgc,
                                               'Wgc': Wgc, 'Lgc': Lgc}
        self.accessible_undetermined_parameters = {'epsilon_gdl': epsilon_gdl, 'tau': tau, 'epsilon_mc': epsilon_mc,
                                                   'epsilon_c': epsilon_c, 'e': e, 'kappa_co': kappa_co, 'Re': Re,
                                                   'i0_c_ref': i0_c_ref, 'kappa_c': kappa_c, 'a_slim': a_slim,
                                                   'b_slim': b_slim, 'a_switch': a_switch, 'C_dl': C_dl}
        self.computing_parameters = {'max_step': max_step, 'n_gdl': n_gdl, 't_purge': t_purge,
                                     'type_fuel_cell': type_fuel_cell, 'type_current': type_current,
                                     'type_auxiliary': type_auxiliary, 'type_control': type_control,
                                     'type_purge': type_purge, 'type_display': type_display, 'type_plot': type_plot}
        self.parameters = {**self.current_parameters, **self.accessible_physical_parameters,
                           **self.accessible_undetermined_parameters, **self.computing_parameters}
        if self.operating_inputs['Pa_des'] < Pext or self.operating_inputs['Pc_des'] < Pext:
            raise ValueError('The desired pressure is too low. It cannot be lower than the pressure outside the stack.')

        # Initialize the variables' dictionary.
        self.solver_variable_names = ['C_v_agc', 'C_v_agdl', 'C_v_acl', 'C_v_ccl', 'C_v_cgdl', 'C_v_cgc', 's_agdl',
                                      's_acl', 's_ccl', 's_cgdl', 'lambda_acl', 'lambda_mem', 'lambda_ccl', 'C_H2_agc',
                                      'C_H2_agdl', 'C_H2_acl', 'C_O2_ccl', 'C_O2_cgdl', 'C_O2_cgc', 'C_N2', 'eta_c',
                                      'Pasm', 'Paem', 'Pcsm', 'Pcem', 'Phi_asm', 'Phi_aem', 'Phi_csm', 'Phi_cem',
                                      'Wcp', 'Wa_inj', 'Wc_inj', 'Abp_a', 'Abp_c']
        self.solver_variable_names_extension()  # Several points are considered in each GDL and must be inserted into
        #                                        the solver_variable_names.
        self.all_variable_names = self.solver_variable_names + ['t', 'Ucell', 'S_sorp_acl', 'S_sorp_ccl'] + \
                                  ['J_lambda_mem_acl', 'J_lambda_mem_ccl', 'Pagc', 'Pcgc', 'Phi_a_des', 'Phi_c_des']
        self.variables = {key: [] for key in self.all_variable_names}

        # Initialize the control_variables dictionary.
        self.control_variables = {'t_control_Phi': self.parameters['t_step'][0],
                                  'Phi_a_des': self.operating_inputs['Phi_a_des'],
                                  'Phi_c_des': self.operating_inputs['Phi_c_des']}

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
                             max_step=self.parameters['max_step'], events=event_negative,
                             args=(self.operating_inputs, self.parameters, self.solver_variable_names,
                                   self.control_variables))

        #       Recover the variable values calculated by the solver into the dictionary.
        self._recovery()

        #       Calculate the cell voltage after computing the internal states of the cell.
        self.variables["Ucell"].extend(calculate_cell_voltage(self.variables, self.operating_inputs, self.parameters))

    def solver_variable_names_extension(self):
        """Several points are considered in each GDL and must be inserted into the solver_variable_names.
        """

        new_points_location = ['C_v_agdl', 'C_v_cgdl', 's_agdl', 's_cgdl', 'C_H2_agdl', 'C_O2_cgdl']
        for variable in new_points_location:
            index = self.solver_variable_names.index(variable)
            # Delete the previous points
            self.solver_variable_names.pop(index)
            # Increase the number of points
            self.solver_variable_names[index:index] = [f'{variable}_{i}' for i in
                                                       range(1, self.parameters['n_gdl'] + 1)]

    def _create_time_interval(self):
        """Calculate the time intervals for numerical resolution, according to the current chosen,
        if it is not provided.

        Returns
        -------
        list
            Time interval for numerical resolution. It is used when initial_variable_values == None.
        """

        # Extraction of the parameters
        t_step, delta_pola, i_max_pola = self.parameters['t_step'], self.parameters['delta_pola'], self.parameters['i_max_pola']
        type_current = self.parameters['type_current']

        # Recovery of the good time interval
        if type_current == "step":
            t0_step, tf_step, delta_t_load_step, delta_t_dyn_step = t_step
            t0_interval = t0_step
            tf_interval = tf_step
        elif type_current == "polarization":
            delta_t_load_pola, delta_t_break_pola, delta_i_pola, delta_t_ini_pola = delta_pola
            t0_interval = 0
            tf_interval = delta_t_ini_pola + int(i_max_pola / delta_i_pola + 1) * (delta_t_load_pola + delta_t_break_pola)
        else:  # EIS time_interval is calculated in the main.py file.
            raise ValueError("Please enter a recognized type_current option for calculating the time interval.")

        # To be reviewed
        self.control_variables['t_control_Phi'] = t0_interval

        return [t0_interval, tf_interval]

    def _create_initial_variable_values(self):
        """Create the initial values of the solver variables if it is not provided.
        It is generated considering an equilibrium at the operating inputs without current.

        Returns
        -------
        list
            Initial values of the solver variables. It is used when initial_variable_values == None.
        """

        # Extraction of the operating inputs and parameters
        current_density, Tfc = self.operating_inputs['current_density'], self.operating_inputs['Tfc']
        Pa_des, Pc_des = self.operating_inputs['Pa_des'], self.operating_inputs['Pc_des']
        Phi_a_des, Phi_c_des = self.operating_inputs['Phi_a_des'], self.operating_inputs['Phi_c_des']
        Hmem, kappa_co, i0_c_ref, = self.parameters['Hmem'], self.parameters['kappa_co'], self.parameters['i0_c_ref']
        kappa_c = self.parameters['kappa_c']
        a_slim, b_slim, a_switch = self.parameters['a_slim'], self.parameters['b_slim'], self.parameters['a_switch']
        n_gdl = self.parameters['n_gdl']

        # Mean value of the operating inputs
        Phi_des_moy = (Phi_a_des + Phi_c_des) / 2
        P_des_moy = (Pa_des + Pc_des) / 2

        # Initial fuel cell states
        #   Intermediate values
        Psat_ini = 101325 * 10 ** (-2.1794 + 0.02953 * (Tfc - 273.15) - 9.1837e-5 * (Tfc - 273.15) ** 2 +
                                   1.4454e-7 * (Tfc - 273.15) ** 3)
        slim = a_slim * (Pc_des / 1e5) + b_slim
        s_switch = a_switch * slim
        #   Initial fuel cell states
        C_v_ini = Phi_des_moy * Psat_ini / (R * Tfc)  # mol.m-3. It is the initial vapor concentration.
        C_H2_ini = (P_des_moy - Phi_des_moy * Psat_ini) / (R * Tfc)  # mol.m-3. It is the initial H2 concentration
        #                                                              in the fuel cell.
        C_O2_ini = yO2_ext * (P_des_moy - Phi_des_moy * Psat_ini) / (R * Tfc)  # mol.m-3. It is the initial O2
        #                                                                        concentration in the fuel cell.
        C_N2_ini = (1 - yO2_ext) * (P_des_moy - Phi_des_moy * Psat_ini) / (R * Tfc)  # mol.m-3. It is the initial N2
        #                                                                              concentration in the fuel cell.
        s_ini = 0  # It is the initial liquid water saturation in the fuel cell.
        lambda_mem_ini = lambda_eq(C_v_ini, s_ini, Tfc, Kshape)  # It is the initial water
        #                                                                                  content in the fuel cell.
        i_fc_ini = current_density(self.time_interval[0], self.parameters)
        i_n_ini = 2 * F * R * Tfc / Hmem * C_H2_ini * k_H2(lambda_mem_ini, Tfc, kappa_co) + \
                  4 * F * R * Tfc / Hmem * C_O2_ini * k_O2(lambda_mem_ini, Tfc, kappa_co)
        f_drop_ini = 0.5 * (1.0 - np.tanh((4 * s_ini - 2 * slim - 2 * s_switch) / (slim - s_switch)))
        eta_c_ini = 1 / f_drop_ini * R * Tfc / (alpha_c * F) * \
                    np.log((i_fc_ini + i_n_ini) / i0_c_ref * (C_O2ref / C_O2_ini) ** kappa_c)  # It is the initial
        #                                                                       cathode overpotential in the fuel cell.

        # Initial auxiliary system state
        Pasm_ini, Paem_ini = Pa_des, Pa_des  # Pa. It is the supply/exhaust manifold pressure at the anode side.
        Pcsm_ini, Pcem_ini = Pc_des, Pc_des  # Pa. It is the supply/exhaust manifold pressure at the cathode side.
        Phi_asm_ini, Phi_aem_ini = Phi_a_des, Phi_a_des  # It is the supply/exhaust manifold relative humidity
        #                                                  at the anode side.
        Phi_csm_ini, Phi_cem_ini = Phi_c_des, Phi_c_des  # It is the supply/exhaust manifold relative humidity
        #                                                  at the cathode side.
        Wcp_ini = 0  # kg.s-1. It is the flow rate of the air compressor.
        Wa_inj_ini = 0  # kg.s-1. It is the flow rate of the air compressor at the anode side.
        Wc_inj_ini = 0  # kg.s-1. It is the flow rate of the air compressor at the cathode side.
        Abp_a_ini = 0  # It is the throttle area of the back pressure valve at the anode.
        Abp_c_ini = 0  # It is the throttle area of the back pressure valve at the cathode.

        # Main variable initialization
        C_v_agc, C_v_agdl, C_v_acl, C_v_ccl, C_v_cgdl, C_v_cgc = [C_v_ini] * 6
        s_agdl, s_acl, s_ccl, s_cgdl = [s_ini] * 4
        s_boundary = 0  # Dirichlet boundary condition
        lambda_acl, lambda_mem, lambda_ccl = [lambda_mem_ini] * 3
        C_H2_agc, C_H2_agdl, C_H2_acl = C_H2_ini, C_H2_ini, C_H2_ini
        C_O2_ccl, C_O2_cgdl, C_O2_cgc = C_O2_ini, C_O2_ini, C_O2_ini
        C_N2, eta_c = C_N2_ini, eta_c_ini
        Pasm, Paem, Pcsm, Pcem = Pasm_ini, Paem_ini, Pcsm_ini, Pcem_ini
        Phi_asm, Phi_aem, Phi_csm, Phi_cem = Phi_asm_ini, Phi_aem_ini, Phi_csm_ini, Phi_cem_ini
        Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c = Wcp_ini, Wa_inj_ini, Wc_inj_ini, Abp_a_ini, Abp_c_ini

        # Gathering of the variables initial value into one list
        initial_variable_values = [C_v_agc] + [C_v_agdl] * n_gdl + [C_v_acl, C_v_ccl] + [C_v_cgdl] * n_gdl + \
                                  [C_v_cgc] + \
                                  [s_boundary] + [s_agdl] * (n_gdl - 1) + [s_acl, s_ccl] + [s_cgdl] * (n_gdl - 1) + \
                                  [s_boundary] + [lambda_acl, lambda_mem, lambda_ccl] + \
                                  [C_H2_agc] + [C_H2_agdl] * n_gdl + [C_H2_acl, C_O2_ccl] + [C_O2_cgdl] * n_gdl + \
                                  [C_O2_cgc, C_N2] + [eta_c] + \
                                  [Pasm, Paem, Pcsm, Pcem, Phi_asm, Phi_aem, Phi_csm, Phi_cem] + \
                                  [Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c]

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
        #   The control variables should be reinitialized. To be reviewed.
        if self.parameters['type_current'] == "step":
            self.control_variables['t_control_Phi'] = self.parameters['t_step'][0]
        else:
            self.control_variables['t_control_Phi'] = 0
        self.control_variables['Phi_a_des'] = self.operating_inputs['Phi_a_des']
        self.control_variables['Phi_c_des'] = self.operating_inputs['Phi_c_des']

        for j in range(len(self.sol.t)):  # For each time...
            # ... recovery of i_fc.
            i_fc = self.operating_inputs["current_density"](self.variables['t'][j], self.parameters)
            # ... recovery of S_sorp_acl, S_sorp_ccl, Jmem_acl, Jmem_ccl, Pagc, Pcgc.
            last_solver_variables = {key: self.variables[key][j] for key in self.solver_variable_names}
            flows_recovery = calculate_flows(self.variables['t'][j], last_solver_variables, self.control_variables,
                                             i_fc, self.operating_inputs, self.parameters)
            for key in ['S_sorp_acl', 'S_sorp_ccl', 'J_lambda_mem_acl', 'J_lambda_mem_ccl', 'Pagc', 'Pcgc']:
                self.variables[key].append(flows_recovery[key])
            # ... recovery of Phi_a_des and Phi_c_des.
            if self.parameters["type_control"] == "Phi_des":
                sv = {'lambda_mem': self.variables['lambda_mem'][j], 's_ccl': self.variables['s_ccl'][j]}
                control_operating_conditions(self.variables['t'][j], sv, self.operating_inputs,
                                             self.parameters, self.control_variables)
                for key in ['Phi_a_des', 'Phi_c_des']: self.variables[key].append(self.control_variables[key])

    def Display(self, ax1=None, ax2=None):
        """Display the plots of the program.

        Parameters
        ----------
        ax1 : matplotlib.axes.Axes, optional
            Axes for the first set of plots. The default is None.
        ax2 : matplotlib.axes.Axes, optional
            Axes for the second set of plots. The default is None.
        """

        # Extraction of the operating inputs and parameters
        Tfc = self.operating_inputs['Tfc']
        n_gdl, type_fuel_cell = self.parameters['n_gdl'], self.parameters['type_fuel_cell']
        type_current, type_display = self.parameters['type_current'], self.parameters['type_display']

        # Parameters' preparation
        n = len(self.variables['t'])
        subfolder_name = type_fuel_cell[:type_fuel_cell.rfind('_')] if type_fuel_cell.rfind('_') != -1 \
            else type_fuel_cell

        # Display
        if type_current == "step":
            if type_display == "multiple":

                figs, axes = zip(*[plt.subplots(figsize=(6, 6)) for _ in range(10)])

                plot_ifc(self.variables, self.operating_inputs, self.parameters, n, axes[0])
                plot_J(self.variables, self.parameters, axes[1])
                plot_C_v(self.variables, n_gdl, C_v_sat(Tfc), n, axes[2])
                plot_lambda(self.variables, self.operating_inputs, self.parameters, axes[3])
                plot_s(self.variables, self.operating_inputs, self.parameters, axes[4])
                plot_C_O2(self.variables, n_gdl, axes[5])
                plot_C_H2(self.variables, n_gdl, axes[6])
                plot_C_N2(self.variables, axes[7])
                plot_Ucell(self.variables, axes[8])
                plot_P(self.variables, axes[9])

                # Considering the number of plots, the saving instructions are made here and not in the main.py file.
                saving_instructions("results", subfolder_name, "step_current_ifc_1.pdf", figs[0])
                saving_instructions("results", subfolder_name, "step_current_J_1.pdf", figs[1])
                saving_instructions("results", subfolder_name, "step_current_Cv_1.pdf", figs[2])
                saving_instructions("results", subfolder_name, "step_current_lambda_1.pdf", figs[3])
                saving_instructions("results", subfolder_name, "step_current_s_1.pdf", figs[4])
                saving_instructions("results", subfolder_name, "step_current_C_O2_1.pdf", figs[5])
                saving_instructions("results", subfolder_name, "step_current_C_H2_1.pdf", figs[6])
                saving_instructions("results", subfolder_name, "step_current_C_N2_1.pdf", figs[7])
                saving_instructions("results", subfolder_name, "step_current_Ucell_1.pdf", figs[8])
                saving_instructions("results", subfolder_name, "step_current_P_1.pdf", figs[9])

                plt.pause(0.001)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "synthetic":

                plot_ifc(self.variables, self.operating_inputs, self.parameters, n, ax1[0, 0])
                plot_Ucell(self.variables, ax1[0, 1])
                plot_J(self.variables, self.parameters, ax1[0, 2])
                plot_C_v(self.variables, n_gdl, C_v_sat(Tfc), n, ax1[1, 0])
                plot_s(self.variables, self.operating_inputs, self.parameters, ax1[1, 1])
                plot_lambda(self.variables, self.operating_inputs, self.parameters, ax1[1, 2])
                plot_C_H2(self.variables, n_gdl, ax1[2, 0])
                plot_C_O2(self.variables, n_gdl, ax1[2, 1])
                plot_P(self.variables, ax1[2, 2])

                plt.pause(0.001)  # A break is necessary to plot the new points in dynamic mode

        elif type_current == "polarization":
            if type_display == "multiple":

                plot_polarisation_curve(self.variables, self.operating_inputs, self.parameters, ax1[0])
                plot_power_density_curve(self.variables, self.operating_inputs, self.parameters, n, ax1[1])
                plot_cell_efficiency(self.variables, self.operating_inputs, self.parameters, n, ax1[2])

                plot_Phi_des(self.variables, self.operating_inputs, self.parameters, ax2[0])
                plot_lambda(self.variables, self.operating_inputs, self.parameters, ax2[1])
                plot_s(self.variables, self.operating_inputs, self.parameters, ax2[2])

                plt.pause(0.001)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "synthetic":

                plot_polarisation_curve(self.variables, self.operating_inputs, self.parameters, ax1)
                plt.pause(0.001)  # A break is necessary to plot the new points in dynamic mode

        elif type_current == "EIS":
            if type_display == "multiple":

                Fourier_results = make_Fourier_transformation(self.variables, self.operating_inputs, self.parameters)
                plot_EIS_curve_Nyquist(self.parameters, Fourier_results, ax1)
                plot_EIS_curve_Bode_amplitude(self.parameters, Fourier_results, ax2[0])
                plot_EIS_curve_Bode_angle(Fourier_results, ax2[1])

                # # Tests to verify the accuracy of EIS simulation.
                # plot_EIS_curve_tests(self.variables, self.operating_inputs, self.parameters, Fourier_results)

                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode

            elif type_display == "synthetic":

                Fourier_results = make_Fourier_transformation(self.variables, self.operating_inputs, self.parameters)
                plot_EIS_curve_Nyquist(self.parameters, Fourier_results, ax1[0])
                plot_EIS_curve_Bode_amplitude(self.parameters, Fourier_results, ax1[1])
                plot_EIS_curve_Bode_angle(Fourier_results, ax1[2])

                # # Tests to verify the accuracy of EIS simulation.
                # plot_EIS_curve_tests(self.variables, self.operating_inputs, self.parameters, Fourier_results)

                plt.pause(0.1)  # A break is necessary to plot the new points in dynamic mode
