# -*- coding: utf-8 -*-

"""This module contains some of the required functions for the GUI.py file.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import time
import numpy as np
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk

# Importing constants' value and functions
from configuration.settings import current_density_parameters, computing_parameters
from model.AlphaPEM import AlphaPEM
from modules.settings_modules import stored_operating_inputs, stored_physical_parameters, EIS_parameters
from modules.main_modules import figures_preparation


# _____________________________________________________GUI modules_____________________________________________________

def changeValue(operating_conditions_frame, accessible_parameters_frame, undetermined_parameters_frame,
                current_density_parameters_frame, computing_parameters_frame, choice_operating_conditions,
                choice_accessible_parameters, choice_undetermined_parameters, choice_current_density_parameters,
                choice_computing_parameters, choices_buttons):
    """This function is called when the user selects a specific option from a dropdown menu for the type of fuel cell.
    Depending on the selected option, it either hides or displays specific input fields (labels or entry widgets) on
    the GUI.

    Parameters
    ----------
    operating_conditions_frame : ttk.Frame
        The frame where the graphical elements for the operating condition and the choice of fuel cell are placed.
    accessible_parameters_frame : ttk.Frame
        The frame where the graphical elements for the accessible physical parameters are placed.
    undetermined_parameters_frame : ttk.Frame
        The frame where the graphical elements for the undetermined physical parameters are placed.
    current_density_parameters_frame : ttk.Frame
        The frame where the graphical elements for the current density parameters are placed.
    computing_parameters_frame : ttk.Frame
        The frame where the graphical elements for the computing parameters are placed.
    choice_operating_conditions : dict
        A dictionary containing the operating condition information.
    choice_accessible_parameters : dict
        A dictionary containing the accessible physical parameter information.
    choice_undetermined_parameters : dict
        A dictionary containing the undetermined physical parameter information.
    choice_current_density_parameters : dict
        A dictionary containing the current density parameter information.
    choice_computing_parameters : dict
        A dictionary containing the computing parameter information.
    choices_buttons : dict
        A dictionary containing the button information.
    """

    if choices_buttons['type_fuel_cell']['value'].get() != 'Enter your specifications':
        # Recovers the new settings
        recover_for_display_operating_inputs_and_physical_parameters(choice_operating_conditions,
                                                                     choice_accessible_parameters,
                                                                     choice_undetermined_parameters,
                                                                     choice_current_density_parameters,
                                                                     choice_computing_parameters, choices_buttons)
        # Display the labels for ...
        #       operating conditions
        for k, v in choice_operating_conditions.items():
            ttk.Label(operating_conditions_frame, width=7, anchor='w', textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       accessible physical parameters
        for k, v in choice_accessible_parameters.items():
            ttk.Label(accessible_parameters_frame, width=7, anchor='w', textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       undetermined physical parameters
        for k, v in choice_undetermined_parameters.items():
            ttk.Label(undetermined_parameters_frame, width=7, anchor='w', textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       current density parameters
        for k, v in choice_current_density_parameters.items():
            ttk.Label(current_density_parameters_frame, width=7, anchor='w', textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       computing parameters
        for k, v in choice_computing_parameters.items():
            ttk.Label(computing_parameters_frame, width=7, anchor='w', textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)

    else:  # choices_buttons['type_fuel_cell']['value'].get() == 'Enter your specifications':
        # Saves and displays the user entries for ...
        #       operating conditions
        for k, v in choice_operating_conditions.items():
            ttk.Entry(operating_conditions_frame, width=7, textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       accessible physical parameters
        for k, v in choice_accessible_parameters.items():
            ttk.Entry(accessible_parameters_frame, width=7, textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       undetermined physical parameters
        for k, v in choice_undetermined_parameters.items():
            ttk.Entry(undetermined_parameters_frame, width=7, textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       current density parameters
        for k, v in choice_current_density_parameters.items():
            ttk.Entry(current_density_parameters_frame, width=7, textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)
        #       computing parameters
        for k, v in choice_computing_parameters.items():
            ttk.Entry(computing_parameters_frame, width=7, textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)


def display_parameter_labels(operating_conditions_frame, accessible_parameters_frame, undetermined_parameters_frame,
                             current_density_parameters_frame, computing_parameters_frame, choice_operating_conditions,
                             choice_accessible_parameters, choice_undetermined_parameters,
                             choice_current_density_parameters, choice_computing_parameters):
    """This function displays labels on the GUI, representing operating conditions and physical parameters, without
    their actual values.

    Parameters
    ----------
    operating_conditions_frame : ttk.Frame
        The frame where the graphical elements for the operating condition and the choice of fuel cell are placed.
    accessible_parameters_frame : ttk.Frame
        The frame where the graphical elements for the accessible physical parameters are placed.
    undetermined_parameters_frame : ttk.Frame
        The frame where the graphical elements for the undetermined physical parameters are placed.
    current_density_parameters_frame : ttk.Frame
        The frame where the graphical elements for the current density parameters are placed.
    computing_parameters_frame : ttk.Frame
        The frame where the graphical elements for the computing parameters are placed.
    choice_operating_conditions : dict
        A dictionary containing the operating condition information.
    choice_accessible_parameters : dict
        A dictionary containing the accessible physical parameter information.
    choice_undetermined_parameters : dict
        A dictionary containing the undetermined physical parameter information.
    choice_current_density_parameters : dict
        A dictionary containing the current density parameter information.
    choice_computing_parameters : dict
        A dictionary containing the computing parameter information.
    """

    # Display the titles
    ttk.Label(operating_conditions_frame, text='Operating conditions', font=('cmr10', 12, 'bold')). \
        grid(row=1, column=0, columnspan=6, ipady=15)
    ttk.Label(accessible_parameters_frame, text='Accessible physical parameters', font=('cmr10', 12, 'bold')). \
        grid(row=0, column=0, columnspan=6, ipady=15)

    # Display the labels for ...
    #       operating conditions
    for k, v in choice_operating_conditions.items():
        ttk.Label(operating_conditions_frame, text=k, font=('cmr10', 10)). \
            grid(row=v['label_row'], column=v['label_column'] - 1, sticky="w")
    #       accessible physical parameters
    for k, v in choice_accessible_parameters.items():
        ttk.Label(accessible_parameters_frame, text=k, font=('cmr10', 10)). \
            grid(row=v['label_row'], column=v['label_column'] - 1, sticky="w")
    #       undetermined physical parameters
    for k, v in choice_undetermined_parameters.items():
        ttk.Label(undetermined_parameters_frame, text=k, font=('cmr10', 10)). \
            grid(row=v['label_row'], column=v['label_column'] - 1, sticky="w")
    #       current density parameters
    ttk.Label(current_density_parameters_frame, text='Step current parameters', font=('cmr10', 10, 'bold')). \
        grid(row=0, column=0, columnspan=2, sticky="w")
    ttk.Label(current_density_parameters_frame, text='Polarization current parameters', font=('cmr10', 10, 'bold')). \
        grid(row=2, column=0, columnspan=2, sticky="w")
    ttk.Label(current_density_parameters_frame, text='EIS current parameters', font=('cmr10', 10, 'bold')). \
        grid(row=5, column=0, columnspan=2, sticky="w")
    for k, v in choice_current_density_parameters.items():
        ttk.Label(current_density_parameters_frame, text=k, font=('cmr10', 10)). \
            grid(row=v['label_row'], column=v['label_column'] - 1, sticky="w")
    #       computing parameters
    for k, v in choice_computing_parameters.items():
        ttk.Label(computing_parameters_frame, text=k, font=('cmr10', 10)). \
            grid(row=v['label_row'], column=v['label_column'] - 1, sticky="w")


def display_parameters_value(operating_conditions_frame, accessible_parameters_frame, undetermined_parameters_frame,
                             current_density_parameters_frame, computing_parameters_frame, choice_operating_conditions,
                             choice_accessible_parameters, choice_undetermined_parameters,
                             choice_current_density_parameters, choice_computing_parameters):
    """This function displays entry widgets on the GUI. There, the user can enter values for operating conditions and
    physical parameters.

    Parameters
    ----------
    operating_conditions_frame : ttk.Frame
        The frame where the graphical elements for the operating condition and the choice of fuel cell are placed.
    accessible_parameters_frame : ttk.Frame
        The frame where the graphical elements for the accessible physical parameters are placed.
    undetermined_parameters_frame : ttk.Frame
        The frame where the graphical elements for the undetermined physical parameters are placed.
    current_density_parameters_frame : ttk.Frame
        The frame where the graphical elements for the current density parameters are placed.
    computing_parameters_frame : ttk.Frame
        The frame where the graphical elements for the computing parameters are placed.
    choice_operating_conditions : dict
        A dictionary containing the operating condition information.
    choice_accessible_parameters : dict
        A dictionary containing the accessible physical parameter information.
    choice_undetermined_parameters : dict
        A dictionary containing the undetermined physical parameter information.
    choice_current_density_parameters : dict
        A dictionary containing the current density parameter information.
    choice_computing_parameters : dict
        A dictionary containing the computing parameter information.
    """
    # Display the value for ...
    #       operating conditions
    for k, v in choice_operating_conditions.items():
        ttk.Entry(operating_conditions_frame, width=7, textvariable=v['value']). \
            grid(row=v['label_row'], column=v['label_column'], padx=5)
    #       accessible physical parameters
    for k, v in choice_accessible_parameters.items():
        ttk.Entry(accessible_parameters_frame, width=7, textvariable=v['value']). \
            grid(row=v['label_row'], column=v['label_column'], padx=5)
    #       undetermined physical parameters
    for k, v in choice_undetermined_parameters.items():
        ttk.Entry(undetermined_parameters_frame, width=7, textvariable=v['value']). \
            grid(row=v['label_row'], column=v['label_column'], padx=5)
    #       current density parameters
    for k, v in choice_current_density_parameters.items():
        ttk.Entry(current_density_parameters_frame, width=7, textvariable=v['value']). \
            grid(row=v['label_row'], column=v['label_column'], padx=5)
    #       computing parameters
    for k, v in choice_computing_parameters.items():
        ttk.Entry(computing_parameters_frame, width=7, textvariable=v['value']). \
            grid(row=v['label_row'], column=v['label_column'], padx=5)


def display_radiobuttons(model_possibilities_frame, choices_buttons):
    """This function displays radiobuttons on the GUI, allowing the user to make choices for control, results display,
    plot style, etc.

    Parameters
    ----------
    model_possibilities_frame : ttk.Frame
        The frame where the graphical elements for the model possibilities and the choice of current density are placed.
    choices_buttons : dict
        A dictionary containing the button information.
    """

    ttk.Label(model_possibilities_frame, text='Model possibilities', font=('cmr10', 12, 'bold')) \
        .grid(row=0, column=0, columnspan=6, ipady=15)

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text='Auxiliaries: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=0, columnspan=1, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='No auxiliaries', value=0,
                    variable=choices_buttons['type_auxiliary']['value']). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Forced-convective cathode\nwith anodic recirculation', value=1,
                    variable=choices_buttons['type_auxiliary']['value']). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=3, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Forced-convective cathode\nwith flow-through anode', value=2,
                    variable=choices_buttons['type_auxiliary']['value']). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=4, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text='Control: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_control']['label_row'], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='No control', value=0,
                    variable=choices_buttons['type_control']['value']). \
        grid(row=choices_buttons['type_control']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Humidity', value=1,
                    variable=choices_buttons['type_control']['value']). \
        grid(row=choices_buttons['type_control']['label_row'], column=3, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text='Purge: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_purge']['label_row'], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='No purge', value=0,
                    variable=choices_buttons['type_purge']['value']). \
        grid(row=choices_buttons['type_purge']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Periodic', value=1,
                    variable=choices_buttons['type_purge']['value']). \
        grid(row=choices_buttons['type_purge']['label_row'], column=3, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Constant', value=2,
                    variable=choices_buttons['type_purge']['value']). \
        grid(row=choices_buttons['type_purge']['label_row'], column=4, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text='Display: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_display']['label_row'], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='No display', value=0,
                    variable=choices_buttons['type_display']['value']). \
        grid(row=choices_buttons['type_display']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Synthetic', value=1,
                    variable=choices_buttons['type_display']['value']). \
        grid(row=choices_buttons['type_display']['label_row'], column=3, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Multiple', value=2,
                    variable=choices_buttons['type_display']['value']). \
        grid(row=choices_buttons['type_display']['label_row'], column=4, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text='Plot: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_plot']['label_row'], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Fixed', value=0,
                    variable=choices_buttons['type_plot']['value']). \
        grid(row=choices_buttons['type_plot']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text='Dynamic', value=1,
                    variable=choices_buttons['type_plot']['value']). \
        grid(row=choices_buttons['type_plot']['label_row'], column=3, sticky="w")


def recover_for_display_operating_inputs_and_physical_parameters(choice_operating_conditions,
                                                                 choice_accessible_parameters,
                                                                 choice_undetermined_parameters,
                                                                 choice_current_density_parameters,
                                                                 choice_computing_parameters, choice_buttons):
    """This function retrieves parameter values for predefined stacks (e.g., "EH-31 1.5 bar (2021)", "Biao Xie 1.0 bar
    (2015)", etc.) and converts them to appropriate units for display on the GUI.

    Parameters
    ----------
    choice_operating_conditions : dict
        A dictionary containing the operating condition information.
    choice_accessible_parameters : dict
        A dictionary containing the accessible physical parameter information.
    choice_undetermined_parameters : dict
        A dictionary containing the undetermined physical parameter information.
    choice_current_density_parameters : dict
        A dictionary containing the current density parameter information.
    choice_computing_parameters : dict
        A dictionary containing the computing parameter information.
    choice_buttons : dict
        A dictionary containing the button information.
    """

    if choice_buttons['type_fuel_cell']['value'].get() == "EH-31 1.5 bar (2021)":
        type_fuel_cell = "EH-31_1.5"
    elif choice_buttons['type_fuel_cell']['value'].get() == "EH-31 2.0 bar (2021)":
        type_fuel_cell = "EH-31_2.0"
    elif choice_buttons['type_fuel_cell']['value'].get() == "EH-31 2.25 bar (2021)":
        type_fuel_cell = "EH-31_2.25"
    elif choice_buttons['type_fuel_cell']['value'].get() == "EH-31 2.5 bar (2021)":
        type_fuel_cell = "EH-31_2.5"
    elif choice_buttons['type_fuel_cell']['value'].get() == "Linhao Fan (2010)":
        type_fuel_cell = "LF"

    (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
     current_density) = current_density_parameters()

    T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, i_max_pola = stored_operating_inputs(type_fuel_cell)

    (Hacl, Hccl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_cl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc,
     Lgc, Vsm, Vem, A_T, Aact, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl) = \
        stored_physical_parameters(type_fuel_cell)

    n_gdl, n_mpl, t_purge, rtol, atol, step_current_parameters = computing_parameters(step_current_parameters, Hgdl, Hmpl, Hacl, type_fuel_cell)

    # operating conditions recovery
    choice_operating_conditions['Temperature - Tfc (°C)']['value'].set(round(T_des - 273.15, 4))  # °C
    choice_operating_conditions['Anode pressure - Pa (bar)']['value'].set(round(Pa_des / 1e5, 4))  # bar
    choice_operating_conditions['Cathode pressure - Pc (bar)']['value'].set(round(Pc_des / 1e5, 4))  # bar
    choice_operating_conditions['Anode stoichiometry - Sa']['value'].set(round(Sa, 4))
    choice_operating_conditions['Cathode stoichiometry - Sc']['value'].set(round(Sc, 4))
    choice_operating_conditions['Anode humidity - Φa']['value'].set(round(Phi_a_des, 4))
    choice_operating_conditions['Cathode humidity - Φc']['value'].set(round(Phi_c_des, 4))
    choice_operating_conditions['Anode inlet H2 ratio - y_H2_in\n(flow-through anode only)']['value'].set(round(y_H2_in, 4))
    # accessible physical parameters recovery
    choice_accessible_parameters['Active area - Aact (cm²)']['value'].set(round(Aact * 1e4, 4))  # cm²
    choice_accessible_parameters['AGC thickness - Hagc (µm)']['value'].set(round(Hagc * 1e6, 4))  # µm
    choice_accessible_parameters['CGC thickness - Hcgc (µm)']['value'].set(round(Hcgc * 1e6, 4))  # µm
    choice_accessible_parameters['AGC width - Wagc (µm)']['value'].set(round(Wagc * 1e6, 4))  # µm
    choice_accessible_parameters['CGC width - Wcgc (µm)']['value'].set(round(Wcgc * 1e6, 4))  # µm
    choice_accessible_parameters['GC cumulated length - Lgc (m)']['value'].set(round(Lgc, 4))  # µm
    choice_accessible_parameters['Supply manifold volume - Vsm (dm³)']['value'].set(round(Vsm * 1e3, 4))  # dm³
    choice_accessible_parameters['Exhaust manifold volume - Vem (dm³)']['value'].set(round(Vem * 1e3, 4))  # dm³
    choice_accessible_parameters['Exhaust manifold throttle area - A_T (cm²)']['value'].set(round(A_T * 1e4, 4))  # cm²
    # undetermined physical parameters recovery
    choice_undetermined_parameters['GDL thickness - Hgdl (µm)']['value'].set(round(Hgdl * 1e6, 4))  # µm
    choice_undetermined_parameters['MPL thickness - Hmpl (µm)']['value'].set(round(Hmpl * 1e6, 4))  # µm
    choice_undetermined_parameters['ACL thickness - Hacl (µm)']['value'].set(round(Hacl * 1e6, 4))  # µm
    choice_undetermined_parameters['CCL thickness - Hccl (µm)']['value'].set(round(Hccl * 1e6, 4))  # µm
    choice_undetermined_parameters['Membrane thickness - Hmem (µm)']['value'].set(round(Hmem * 1e6, 4))  # µm
    choice_undetermined_parameters['GDL porosity - ε_gdl']['value'].set(round(epsilon_gdl, 4))
    choice_undetermined_parameters['CL porosity - ε_cl']['value'].set(round(epsilon_cl, 4))
    choice_undetermined_parameters['MPL porosity - ε_mpl']['value'].set(round(epsilon_mpl, 4))
    choice_undetermined_parameters['Ionomer volume fraction - ε_mc']['value'].set(round(epsilon_mc, 4))
    choice_undetermined_parameters['Compression ratio - ε_c']['value'].set(round(epsilon_c, 4))
    choice_undetermined_parameters['Capillary exponent - e']['value'].set(e)
    choice_undetermined_parameters['Reference exchange current\ndensity - i0_c_ref (A/m²)']['value'].set(round(i0_c_ref, 4))  # A.m-2
    choice_undetermined_parameters['Crossover correction coefficient\n- κ_co (mol/(m.s.Pa))']['value'].set(round(kappa_co, 4))  # mol.m-1.s-1.Pa-1
    choice_undetermined_parameters['Overpotential correction\nexponent - κ_c']['value'].set(round(kappa_c, 4))
    choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_slim']['value'].set(round(a_slim, 7))
    choice_undetermined_parameters['Limit liquid saturation\ncoefficient - b_slim']['value'].set(round(b_slim, 7))
    choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_switch']['value'].set(round(a_switch, 7))
    choice_undetermined_parameters['Volumetric space-charge layer\ncapacitance - C_scl (F/cm³)']['value'].set(round(C_scl * 1e-6, 4))  # F.cm-3
    # i_max_pola recovery
    choice_current_density_parameters['Maximum current density\n- i_max_pola (A/cm²)']['value'].set(round(i_max_pola / 1e4, 4))  # A/cm²
    # computing parameters recovery
    choice_computing_parameters['Number of GDL nodes - n_gdl']['value'].set(n_gdl)
    choice_computing_parameters['Number of MPL nodes - n_mpl']['value'].set(n_mpl)
    choice_computing_parameters['Solver relative tolerance - rtol']['value'].set(rtol)
    choice_computing_parameters['Solver absolute tolerance - atol']['value'].set(atol)


def recover_for_use_operating_inputs_and_physical_parameters(choice_operating_conditions, choice_accessible_parameters,
                                                             choice_undetermined_parameters,
                                                             choice_current_density_parameters,
                                                             choice_computing_parameters, choice_buttons):
    """This function retrieves and converts the parameter values from the GUI into standard units for further
    calculations.

    Parameters
    ----------
    choice_operating_conditions : dict
        A dictionary containing the operating condition information.
    choice_accessible_parameters : dict
        A dictionary containing the accessible physical parameter information.
    choice_undetermined_parameters : dict
        A dictionary containing the undetermined physical parameter information.
    choice_current_density_parameters : dict
        A dictionary containing the current density parameter information.
    choice_computing_parameters : dict
        A dictionary containing the computing parameter information.
    choice_buttons : dict
        A dictionary containing the button information.
    """
    # operating conditions
    T_des = choice_operating_conditions['Temperature - Tfc (°C)']['value'].get() + 273.15  # K
    Pa_des = choice_operating_conditions['Anode pressure - Pa (bar)']['value'].get() * 1e5  # Pa
    Pc_des = choice_operating_conditions['Cathode pressure - Pc (bar)']['value'].get() * 1e5  # Pa
    Sa = choice_operating_conditions['Anode stoichiometry - Sa']['value'].get()
    Sc = choice_operating_conditions['Cathode stoichiometry - Sc']['value'].get()
    Phi_a_des = choice_operating_conditions['Anode humidity - Φa']['value'].get()
    Phi_c_des = choice_operating_conditions['Cathode humidity - Φc']['value'].get()
    y_H2_in = choice_operating_conditions['Anode inlet H2 ratio - y_H2_in\n(flow-through anode only)']['value'].get()
    # accessible physical parameters
    Aact = choice_accessible_parameters['Active area - Aact (cm²)']['value'].get() * 1e-4  # m²
    Hagc = choice_accessible_parameters['AGC thickness - Hagc (µm)']['value'].get() * 1e-6  # m
    Hcgc = choice_accessible_parameters['CGC thickness - Hcgc (µm)']['value'].get() * 1e-6  # m
    Wagc = choice_accessible_parameters['AGC width - Wagc (µm)']['value'].get() * 1e-6  # m
    Wcgc = choice_accessible_parameters['CGC width - Wcgc (µm)']['value'].get() * 1e-6  # m
    Lgc = choice_accessible_parameters['GC cumulated length - Lgc (m)']['value'].get()  # m
    Vsm = choice_accessible_parameters['Supply manifold volume - Vsm (dm³)']['value'].get() * 1e-3  # m³
    Vem = choice_accessible_parameters['Exhaust manifold volume - Vem (dm³)']['value'].get() * 1e-3 # m³
    A_T = choice_accessible_parameters['Exhaust manifold throttle area - A_T (cm²)']['value'].get() * 1e-4  # m²
    # undetermined physical parameters
    Hgdl = choice_undetermined_parameters['GDL thickness - Hgdl (µm)']['value'].get() * 1e-6  # m
    Hmpl = choice_undetermined_parameters['MPL thickness - Hmpl (µm)']['value'].get() * 1e-6  # m
    Hacl = choice_undetermined_parameters['ACL thickness - Hacl (µm)']['value'].get() * 1e-6  # m
    Hccl = choice_undetermined_parameters['CCL thickness - Hccl (µm)']['value'].get() * 1e-6  # m
    Hmem = choice_undetermined_parameters['Membrane thickness - Hmem (µm)']['value'].get() * 1e-6  # m
    epsilon_gdl = choice_undetermined_parameters['GDL porosity - ε_gdl']['value'].get()
    epsilon_cl = choice_undetermined_parameters['CL porosity - ε_cl']['value'].get()
    epsilon_mpl = choice_undetermined_parameters['MPL porosity - ε_mpl']['value'].get()
    epsilon_mc = choice_undetermined_parameters['Ionomer volume fraction - ε_mc']['value'].get()
    epsilon_c = choice_undetermined_parameters['Compression ratio - ε_c']['value'].get()
    e = choice_undetermined_parameters['Capillary exponent - e']['value'].get()
    i0_c_ref = choice_undetermined_parameters['Reference exchange current\ndensity - i0_c_ref (A/m²)']['value'].get()  # A.m-2
    kappa_co = choice_undetermined_parameters['Crossover correction coefficient\n- κ_co (mol/(m.s.Pa))']['value'].get()  # mol.m-1.s-1.Pa-1
    kappa_c = choice_undetermined_parameters['Overpotential correction\nexponent - κ_c']['value'].get()
    a_slim = choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_slim']['value'].get()
    b_slim = choice_undetermined_parameters['Limit liquid saturation\ncoefficient - b_slim']['value'].get()
    a_switch = choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_switch']['value'].get()
    C_scl = choice_undetermined_parameters['Volumetric space-charge layer\ncapacitance - C_scl (F/cm³)']['value'].get() * 1e6  # F.m-3
    # current density parameters
    delta_t_ini_step = choice_current_density_parameters['Stabilisation time\n- Δt_ini_step (min)']['value'].get() * 60 #s
    delta_t_load_step = choice_current_density_parameters['Loading time\n- Δt_load_step (s)']['value'].get() #s
    delta_t_break_step = choice_current_density_parameters['Breaking time\n- Δt_break_step (min)']['value'].get() * 60 #s
    i_step = choice_current_density_parameters['Current density step\n- i_step (A/cm²)']['value'].get() * 1e4 # A.m-2
    delta_t_dyn_step = choice_computing_parameters['Time for dynamic\ndisplay - Δt_dyn_step (s)']['value'].get() #s
    step_current_parameters = {'delta_t_ini_step': delta_t_ini_step, 'delta_t_load_step': delta_t_load_step,
                               'delta_t_break_step': delta_t_break_step, 'i_step': i_step,
                               'delta_t_dyn_step': delta_t_dyn_step}
    delta_t_ini_pola = choice_current_density_parameters['Stabilisation time\n- Δt_ini_pola (min)']['value'].get() * 60 #s
    delta_t_load_pola = choice_current_density_parameters['Loading time\n- Δt_load_pola (s)']['value'].get() #s
    delta_t_break_pola = choice_current_density_parameters['Breaking time\n- Δt_break_pola (min)']['value'].get() * 60 #s
    delta_i_pola = choice_current_density_parameters['Current density step\n- Δi_pola (A/cm²)']['value'].get() * 1e4 # A.m-2
    i_max_pola = choice_current_density_parameters['Maximum current density\n- i_max_pola (A/cm²)']['value'].get() * 1e4 # A.m-2
    pola_current_parameters = {'delta_t_ini_pola': delta_t_ini_pola, 'delta_t_load_pola': delta_t_load_pola,
                               'delta_t_break_pola': delta_t_break_pola, 'delta_i_pola': delta_i_pola,
                               'i_max_pola': i_max_pola}
    pola_current_for_cali_parameters = None # Calibration is not implemented in the GUI.
    i_EIS = choice_current_density_parameters['Static current\n- i_EIS (A/cm²)']['value'].get() * 1e4  # (A.m-2)
    ratio_EIS = choice_current_density_parameters['Current ratio\n- ratio_EIS (%)']['value'].get() / 100
    f_EIS = (choice_current_density_parameters['Power of the\ninitial frequency\n- f_power_min_EIS']['value'].get(),
             choice_current_density_parameters['Power of the\nfinal frequency\n- f_power_max_EIS']['value'].get(),
             choice_current_density_parameters['Number of frequencies\ntested - nb_f_EIS']['value'].get(),
             choice_current_density_parameters['Number of points\ncalculated - nb_points_EIS']['value'].get())
    t_EIS = EIS_parameters(f_EIS)  # Time parameters for the EIS_current density function.
    # computing parameters
    t_purge = choice_computing_parameters['Purge time - t_purge (s)']['value'].get()  # s
    delta_t_purge = choice_computing_parameters['Time between two purges\n- Δt_purge (s)']['value'].get()  # s
    n_gdl = choice_computing_parameters['Number of GDL nodes - n_gdl']['value'].get()
    n_mpl = choice_computing_parameters['Number of MPL nodes - n_mpl']['value'].get()
    rtol = choice_computing_parameters['Solver relative tolerance - rtol']['value'].get()
    atol = choice_computing_parameters['Solver absolute tolerance - atol']['value'].get()

    if choice_buttons['type_fuel_cell']['value'].get() == "EH-31 1.5 bar (2021)":
        type_fuel_cell = "EH-31_1.5"
    elif choice_buttons['type_fuel_cell']['value'].get() == "EH-31 2.0 bar (2021)":
        type_fuel_cell = "EH-31_2.0"
    elif choice_buttons['type_fuel_cell']['value'].get() == "EH-31 2.25 bar (2021)":
        type_fuel_cell = "EH-31_2.25"
    elif choice_buttons['type_fuel_cell']['value'].get() == "EH-31 2.5 bar (2021)":
        type_fuel_cell = "EH-31_2.5"
    elif choice_buttons['type_fuel_cell']['value'].get() == "Linhao Fan (2010)":
        type_fuel_cell = "LF"
    elif choice_buttons['type_fuel_cell']['value'].get() == "Enter your specifications":
        type_fuel_cell = "manual_setup"

    if choice_buttons['type_auxiliary']['value'].get() == 0:
        type_auxiliary = "no_auxiliary"
    elif choice_buttons['type_auxiliary']['value'].get() == 1:
        type_auxiliary = "forced-convective_cathode_with_anodic_recirculation"
    else:
        type_auxiliary = "forced-convective_cathode_with_flow-through_anode"

    if choice_buttons['type_control']['value'].get() == 0:
        type_control = "no_control"
    else:
        type_control = "Phi_des"

    if choice_buttons['type_purge']['value'].get() == 0:
        type_purge = "no_purge"
    elif choice_buttons['type_purge']['value'].get() == 1:
        type_purge = "periodic_purge"
    else:
        type_purge = "constant_purge"

    if choice_buttons['type_display']['value'].get() == 0:
        type_display = "no_display"
    elif choice_buttons['type_display']['value'].get() == 1:
        type_display = "synthetic"
    else:
        type_display = "multiple"

    if choice_buttons['type_plot']['value'].get() == 0:
        type_plot = "fixed"
    else:
        type_plot = "dynamic"

    return (T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, Aact, Hgdl, Hmpl, Hacl, Hccl, Hmem, Hagc, Hcgc, Wagc,
            Wcgc, Lgc, epsilon_gdl, epsilon_cl, epsilon_mpl, epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c,
            a_slim, b_slim, a_switch, C_scl, step_current_parameters, pola_current_parameters,
            pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS, t_purge, delta_t_purge, n_gdl, n_mpl,
            rtol, atol, type_fuel_cell, type_auxiliary, type_control, type_purge, type_display, type_plot)


def value_control(choice_operating_conditions, choice_accessible_parameters, choice_undetermined_parameters,
                  choice_current_density_parameters, choice_computing_parameters, choice_buttons, current_button):
    """This function checks the integrity of the values entered by the user and returns an empty tuple if they are not
    valid.

    Parameters
    ----------
    choice_operating_conditions : dict
        A dictionary containing the operating condition information.
    choice_accessible_parameters : dict
        A dictionary containing the accessible physical parameter information.
    choice_undetermined_parameters : dict
        A dictionary containing the undetermined physical parameter information.
    choice_current_density_parameters : dict
        A dictionary containing the current density parameter information.
    choice_computing_parameters : dict
        A dictionary containing the computing parameter information.
    choice_buttons : dict
        A dictionary containing the button information.
    current_button : dict
        A dictionary representing the clicked button.
    """

    # The values entered by the user are checked for compliance
    if choice_operating_conditions['Temperature - Tfc (°C)']['value'].get() < 0:
        messagebox.showerror(title='Temperatures', message='Negative temperatures do not exist in the Kelvin scale.')
        choices.clear()
        return
    if choice_operating_conditions['Anode pressure - Pa (bar)']['value'].get() < 0 or \
            choice_operating_conditions['Cathode pressure - Pc (bar)']['value'].get() < 0 or \
            choice_operating_conditions['Cathode pressure - Pc (bar)']['value'].get() > 5.0 or \
            choice_operating_conditions['Cathode pressure - Pc (bar)']['value'].get() > 5.0:
        messagebox.showerror(title='Desired pressures', message='Desired pressure should be positive and bellow 5.0 '
                                                                'bars.')
        choices.clear()
        return
    if choice_operating_conditions['Anode stoichiometry - Sa']['value'].get() < 1 or \
            choice_operating_conditions['Anode stoichiometry - Sa']['value'].get() > 5 or \
            choice_operating_conditions['Cathode stoichiometry - Sc']['value'].get() < 1 or \
            choice_operating_conditions['Cathode stoichiometry - Sc']['value'].get() > 5:
        messagebox.showerror(title='Stoichiometric ratios', message='The stoichiometric ratios Sa and Sc should be '
                                                                    'between 1 and 5.')
        choices.clear()
        return
    if choice_operating_conditions['Anode humidity - Φa']['value'].get() < 0 or \
            choice_operating_conditions['Anode humidity - Φa']['value'].get() > 1 or \
            choice_operating_conditions['Cathode humidity - Φc']['value'].get() < 0 or \
            choice_operating_conditions['Cathode humidity - Φc']['value'].get() > 1:
        messagebox.showerror(title='Desired humidity', message='The desired humidities should be between 0 and 1.')
        choices.clear()
        return
    if choice_operating_conditions['Anode inlet H2 ratio - y_H2_in\n(flow-through anode only)']['value'].get() < 0 or \
            choice_operating_conditions['Anode inlet H2 ratio - y_H2_in\n(flow-through anode only)']['value'].get() > 1:
        messagebox.showerror(title='Anode inlet H2 ratio', message='The anode inlet H2 ratio should be between 0 and 1.')
        choices.clear()
        return
    if choice_accessible_parameters['Active area - Aact (cm²)']['value'].get() < 0:
        messagebox.showerror(title='Active area', message='Negative active area is impossible.')
        choices.clear()
        return
    if choice_accessible_parameters['Supply manifold volume - Vsm (dm³)']['value'].get() < 0 or \
            choice_accessible_parameters['Exhaust manifold volume - Vem (dm³)']['value'].get() < 0 or \
            choice_accessible_parameters['Exhaust manifold throttle area - A_T (cm²)']['value'].get() < 0:
        messagebox.showerror(title='Manifold parameters', message='Negative volumes or area are impossible.')
        choices.clear()
        return
    if choice_accessible_parameters['AGC thickness - Hagc (µm)']['value'].get() < 10 or \
            choice_accessible_parameters['AGC thickness - Hagc (µm)']['value'].get() > 10000 or \
            choice_accessible_parameters['CGC thickness - Hcgc (µm)']['value'].get() < 10 or \
            choice_accessible_parameters['CGC thickness - Hcgc (µm)']['value'].get() > 10000 or \
            choice_accessible_parameters['AGC width - Wagc (µm)']['value'].get() < 10 or \
            choice_accessible_parameters['AGC width - Wagc (µm)']['value'].get() > 10000 or \
            choice_accessible_parameters['CGC width - Wcgc (µm)']['value'].get() < 10 or \
            choice_accessible_parameters['CGC width - Wcgc (µm)']['value'].get() > 10000 or \
            choice_accessible_parameters['GC cumulated length - Lgc (m)']['value'].get() < 0 or \
            choice_accessible_parameters['GC cumulated length - Lgc (m)']['value'].get() > 100:
        messagebox.showerror(title='GC distances', message='GC generally have a thickness and a width between 10µm and '
                                                           '10mm. Also, GC length is generally between 0 and 100m')
        choices.clear()
        return
    if choice_undetermined_parameters['GDL thickness - Hgdl (µm)']['value'].get() < 1 or \
            choice_undetermined_parameters['GDL thickness - Hgdl (µm)']['value'].get() > 1000 or \
            choice_undetermined_parameters['MPL thickness - Hmpl (µm)']['value'].get() < 1 or \
            choice_undetermined_parameters['MPL thickness - Hmpl (µm)']['value'].get() > 1000 or \
            choice_undetermined_parameters['ACL thickness - Hacl (µm)']['value'].get() < 1 or \
            choice_undetermined_parameters['ACL thickness - Hacl (µm)']['value'].get() > 1000 or \
            choice_undetermined_parameters['CCL thickness - Hccl (µm)']['value'].get() < 1 or \
            choice_undetermined_parameters['CCL thickness - Hccl (µm)']['value'].get() > 1000 or \
            choice_undetermined_parameters['Membrane thickness - Hmem (µm)']['value'].get() < 1 or \
            choice_undetermined_parameters['Membrane thickness - Hmem (µm)']['value'].get() > 1000:
        messagebox.showerror(title='MEA thickness', message='All MEA components generally have a thickness between '
                                                            '1µm and 1mm.')
        choices.clear()
        return
    if choice_undetermined_parameters['GDL porosity - ε_gdl']['value'].get() < 0.50 or \
            choice_undetermined_parameters['GDL porosity - ε_gdl']['value'].get() > 0.90:
        messagebox.showerror(title='GDL porosity', message='GDL porosity should be between 0.50 and 0.90.')
        choices.clear()
        return
    if choice_undetermined_parameters['CL porosity - ε_cl']['value'].get() < 0.12 or \
            choice_undetermined_parameters['CL porosity - ε_cl']['value'].get() > 0.60:
        messagebox.showerror(title='GDL porosity', message='CL porosity should be between 0.12 and 0.60.')
        choices.clear()
        return
    if choice_undetermined_parameters['MPL porosity - ε_mpl']['value'].get() < 0.30 or \
            choice_undetermined_parameters['MPL porosity - ε_mpl']['value'].get() > 0.70:
        messagebox.showerror(title='MPL porosity', message='MPL porosity should be between 0.30 and 0.70.')
        choices.clear()
        return
    if choice_undetermined_parameters['Ionomer volume fraction - ε_mc']['value'].get() < 0 or \
            choice_undetermined_parameters['Ionomer volume fraction - ε_mc']['value'].get() > 1:
        messagebox.showerror(title='Ionomer volume fraction', message='Ionomer volume fraction should be between 0 and 1.')
        choices.clear()
        return
    if choice_undetermined_parameters['Compression ratio - ε_c']['value'].get() < 0 or \
            choice_undetermined_parameters['Compression ratio - ε_c']['value'].get() > 1:
        messagebox.showerror(title='Compression ratio', message='The compression ratio should be between 0 and 1.')
        choices.clear()
        return
    if choice_undetermined_parameters['Capillary exponent - e']['value'].get() < 3 or choice_undetermined_parameters['Capillary exponent - e']['value'].get() > 5:
        messagebox.showerror(title='Capillary exponent', message='The capillary exponent should be between 3 and 5 and '
                                                                 'being an integer.')
        choices.clear()
        return
    if choice_undetermined_parameters['Reference exchange current\ndensity - i0_c_ref (A/m²)']['value'].get() < 0.001 or \
            choice_undetermined_parameters['Reference exchange current\ndensity - i0_c_ref (A/m²)']['value'].get() > 500:
        messagebox.showerror(title='Referenced exchange current density', message='The referenced exchange current '
                                                                                  'density is generally between 0.001 '
                                                                                  'and 500 A.m-2.')
        choices.clear()
        return
    if choice_undetermined_parameters['Crossover correction coefficient\n- κ_co (mol/(m.s.Pa))']['value'].get() < 0.01 or \
            choice_undetermined_parameters['Crossover correction coefficient\n- κ_co (mol/(m.s.Pa))']['value'].get() > 100:
        messagebox.showerror(title='Crossover correction coefficient', message='The crossover correction coefficient is'
                                                                               ' generally between 0.01 and 100 '
                                                                               'mol.m-1.s-1.Pa-1.')
        choices.clear()
        return
    if choice_undetermined_parameters['Overpotential correction\nexponent - κ_c']['value'].get() < 0 or \
            choice_undetermined_parameters['Overpotential correction\nexponent - κ_c']['value'].get() > 100:
        messagebox.showerror(title='Overpotential correction exponent', message='The overpotential correction exponent '
                                                                                'is generally between 0 and 100.')
        choices.clear()
        return
    if choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_slim']['value'].get() < 0 or \
            choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_slim']['value'].get() > 1:
        messagebox.showerror(title='Slop of slim function', message='The slop of slim function is generally between 0 '
                                                                    'and 1.')
        choices.clear()
        return
    if choice_undetermined_parameters['Limit liquid saturation\ncoefficient - b_slim']['value'].get() < 0 or \
            choice_undetermined_parameters['Limit liquid saturation\ncoefficient - b_slim']['value'].get() > 1:
        messagebox.showerror(title='Intercept of slim function', message='The intercept of slim function is generally '
                                                                         'between 0 and 1.')
        choices.clear()
        return
    if choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_switch']['value'].get() < 0 or \
            choice_undetermined_parameters['Limit liquid saturation\ncoefficient - a_switch']['value'].get() > 1:
        messagebox.showerror(title='Slop of switch function', message='The slop of switch function is generally between'
                                                                      ' 0 and 1.')
        choices.clear()
        return
    if choice_undetermined_parameters['Volumetric space-charge layer\ncapacitance - C_scl (F/cm³)']['value'].get() < 5 or \
            choice_undetermined_parameters['Volumetric space-charge layer\ncapacitance - C_scl (F/cm³)']['value'].get() > 100:
        messagebox.showerror(title='Double layer capacitance', message='I have not settled yet a range for C_scl.')
        choices.clear()
        return
    if choice_current_density_parameters['Stabilisation time\n- Δt_ini_step (min)']['value'].get() < 0 or \
            choice_current_density_parameters['Loading time\n- Δt_load_step (s)']['value'].get() < 0 or \
            choice_current_density_parameters['Breaking time\n- Δt_break_step (min)']['value'].get() < 0 or \
            choice_computing_parameters['Time for dynamic\ndisplay - Δt_dyn_step (s)']['value'].get() < 0 or \
            choice_current_density_parameters['Stabilisation time\n- Δt_ini_pola (min)']['value'].get() < 0 or \
            choice_current_density_parameters['Loading time\n- Δt_load_pola (s)']['value'].get() < 0 or \
            choice_current_density_parameters['Breaking time\n- Δt_break_pola (min)']['value'].get() < 0 :
        messagebox.showerror(title='Times', message='The times should be positive, t0_step < tf_step and '
                                                    'delta_t_load_step < (tf_step - t0_step).')
        choices.clear()
        return
    if choice_current_density_parameters['Maximum current density\n- i_max_pola (A/cm²)']['value'].get() < 0 or \
            choice_current_density_parameters['Current density step\n- Δi_pola (A/cm²)']['value'].get() < 0 or \
            choice_current_density_parameters['Static current\n- i_EIS (A/cm²)']['value'].get() < 0 or \
            choice_current_density_parameters['Current density step\n- Δi_pola (A/cm²)']['value'].get() > \
            choice_current_density_parameters['Maximum current density\n- i_max_pola (A/cm²)']['value'].get():
        messagebox.showerror(title='Current densities', message='The current densities should be positive, '
                                                                'delta_i_pola < i_max_pola and '
                                                                'i_ini_step < i_final_step.')
        choices.clear()
        return
    if choice_current_density_parameters['Current ratio\n- ratio_EIS (%)']['value'].get() < 0 or \
            choice_current_density_parameters['Current ratio\n- ratio_EIS (%)']['value'].get() > 20:
        messagebox.showerror(title='Ratio EIS', message='Ratio EIS is a percentage of i_EIS and should be between 0 '
                                                        'and 20 for plotting correct EIS.')
        choices.clear()
        return

    if choice_current_density_parameters['Number of frequencies\ntested - nb_f_EIS']['value'].get() < 0 or \
            choice_current_density_parameters['Number of points\ncalculated - nb_points_EIS']['value'].get() < 0 or \
            type(choice_current_density_parameters['Power of the\ninitial frequency\n- f_power_min_EIS']['value'].get()) != int or \
            type(choice_current_density_parameters['Power of the\nfinal frequency\n- f_power_max_EIS']['value'].get()) != int or \
            type(choice_current_density_parameters['Number of frequencies\ntested - nb_f_EIS']['value'].get()) != int or \
            type(choice_current_density_parameters['Number of points\ncalculated - nb_points_EIS']['value'].get()) != int:
        messagebox.showerror(title='f EIS', message='f_EIS parameters should be integer and number of points should '
                                                    'be positive.')
        choices.clear()
        return

    if choice_computing_parameters['Purge time - t_purge (s)']['value'].get() < 0 or \
            choice_computing_parameters['Time between two purges\n- Δt_purge (s)']['value'].get() < 0:
        messagebox.showerror(title='Purge times', message='Negative times does not characterise purges.')
        choices.clear()
        return

    if choice_computing_parameters['Number of GDL nodes - n_gdl']['value'].get() < 1 or \
            type(choice_computing_parameters['Number of GDL nodes - n_gdl']['value'].get()) != int:
        messagebox.showerror(title='n gdl', message='The n_gdl value should be an integer bigger or equal to 1.')
        choices.clear()
        return

    if choice_computing_parameters['Number of MPL nodes - n_mpl']['value'].get() < 1 or \
            type(choice_computing_parameters['Number of MPL nodes - n_mpl']['value'].get()) != int:
        messagebox.showerror(title='n gdl', message='The n_mpl value should be an integer bigger or equal to 1.')
        choices.clear()
        return

    if choice_computing_parameters['Solver relative tolerance - rtol']['value'].get() > 1e-3 or \
            choice_computing_parameters['Solver absolute tolerance - atol']['value'].get() > 1e-3:
        messagebox.showerror(title='Solver tolerance', message='rtol and atol should be lower than 1e-3 to limit the'
                                                               ' numerical errors.')
        choices.clear()
        return

    if current_button == 0 and choice_buttons['type_display']['value'].get() == 2 \
                           and choice_buttons['type_plot']['value'].get() == 1 :
        messagebox.showerror(title='n gdl', message='dynamic plot is not thought to be used with step current and '
                                                    'multiple display. There would be too much plots to handle.')
        choices.clear()
        return


def set_equal_width(frame1, frame2, frame3, frame4, frame5, frame6):
    """
    Adjusts the width of the frames to be equal based on their maximum width.

    Parameters
    ----------
    frame1 : ttk.Frame
        The first frame to be resized.
    frame2 : ttk.Frame
        The second frame to be resized.
    frame3 : ttk.Frame
        The third frame to be resized.
    frame4 : ttk.Frame
        The fourth frame to be resized.
    frame5 : ttk.Frame
        The fifth frame to be resized.
    frame6 : ttk.Frame
        The sixth frame to be resized.
    """

    # Initialisation of the list of widths
    widths = []

    for frame in [frame1, frame2, frame3, frame4, frame5, frame6]:
        # Update the frame sizes
        frame.update_idletasks()
        # Get the current width of all frames
        widths.append(frame.winfo_width())

    # Set all frames to the maximum width
    for frame in [frame1, frame2, frame3, frame4, frame5, frame6]:
        for i in range(6):
            frame.grid_columnconfigure(i, minsize=max(widths) / 5.5)  # Set minimum width of all column to max_width / 5


def launch_AlphaPEM_for_step_current(operating_inputs, current_parameters, accessible_physical_parameters,
                                     undetermined_physical_parameters, computing_parameters):
    """Launch the AlphaPEM simulator for a step current density and display the results.

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
            - epsilon_cl : float
                Anode/cathode CL porosity (undetermined physical parameter).
            - epsilon_mc : float
                Volume fraction of ionomer in the CL (undetermined physical parameter).
            - epsilon_c : float
                Compression ratio of the GDL (undetermined physical parameter).
            - e : float
                Capillary exponent (undetermined physical parameter).
            - i0_c_ref : float
                Reference exchange current density at the cathode in A.m-2 (undetermined physical parameter).
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
            - n_gdl : int
                Number of points considered in the GDL (computing parameter).
            - n_mpl : int
                Number of points considered in the MPL (computing parameter).
            - t_purge : tuple
                Time parameters for purging the system (computing parameter).
                It is the purge time interval 'purge_time' and the time between two purges 'delta_purge'.
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
    """

    # Starting time
    start_time = time.time()

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if computing_parameters['type_plot'] == "dynamic":
        # Check if the type_fuel_cell and type_current are valid
        if computing_parameters['type_current'] == "step" and computing_parameters['type_display'] == "multiple":
            raise ValueError('dynamic plot is not thought to be used with step current and multiple display.' +
                             'There would be too much plots to handle.')

        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = None
        #           Extraction of the parameters
        tf_step = (current_parameters['step_current_parameters']['delta_t_ini_step'] +
                   current_parameters['step_current_parameters']['delta_t_load_step'] +
                   current_parameters['step_current_parameters']['delta_t_break_step'])  # (s).
        delta_t_dyn_step = current_parameters['step_current_parameters']['delta_t_dyn_step']  # (s).
        #           Calculation
        n = int(tf_step / delta_t_dyn_step)  # It is the plot update number.
        time_interval = [0, delta_t_dyn_step]  # (s). It is the initial time interval.

        # Dynamic simulation
        for i in range(n):
            Simulator = AlphaPEM(operating_inputs, current_parameters, accessible_physical_parameters,
                                 undetermined_physical_parameters, computing_parameters, initial_variable_values,
                                 time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                t0_interval = Simulator.variables['t'][-1]
                tf_interval = (i + 2) * delta_t_dyn_step
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator.solver_variable_names:
                initial_variable_values.append(Simulator.variables[x][-1])

            # Display
            if computing_parameters['type_display'] != "no_display":
                Simulator.Display(ax1, ax2, ax3)

    else:  # elif computing_parameters['type_plot'] == "fixed":
        # Simulation
        Simulator = AlphaPEM(operating_inputs, current_parameters, accessible_physical_parameters,
                             undetermined_physical_parameters, computing_parameters)
        # Display
        if computing_parameters['type_display'] != "no_display":
            Simulator.Display(ax1, ax2, ax3)

    # Plot saving
    Simulator.Save_plot(fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)


def launch_AlphaPEM_for_polarization_current(operating_inputs, current_parameters, accessible_physical_parameters,
                                             undetermined_physical_parameters, computing_parameters):
    """Launch the AlphaPEM simulator for a polarization current density and display the results.

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
            - epsilon_cl : float
                Anode/cathode CL porosity (undetermined physical parameter).
            - epsilon_mc : float
                Volume fraction of ionomer in the CL (undetermined physical parameter).
            - epsilon_c : float
                Compression ratio of the GDL (undetermined physical parameter).
            - e : float
                Capillary exponent (undetermined physical parameter).
            - i0_c_ref : float
                Reference exchange current density at the cathode in A.m-2 (undetermined physical parameter).
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
            - n_gdl : int
                Number of points considered in the GDL (computing parameter).
            - n_mpl : int
                Number of points considered in the MPL (computing parameter).
            - t_purge : tuple
                Time parameters for purging the system (computing parameter).
                It is the purge time interval 'purge_time' and the time between two purges 'delta_purge'.
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
    """

    # Starting time
    start_time = time.time()

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if computing_parameters['type_plot'] == "dynamic":
        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = None
        #           Extraction of the parameters
        delta_t_ini_pola = current_parameters['pola_current_parameters']['delta_t_ini_pola']  # (s).
        delta_t_load_pola = current_parameters['pola_current_parameters']['delta_t_load_pola']  # (s).
        delta_t_break_pola = current_parameters['pola_current_parameters']['delta_t_break_pola']  # (s).
        delta_i_pola = current_parameters['pola_current_parameters']['delta_i_pola']  # (A.m-2).
        i_max_pola = current_parameters['pola_current_parameters']['i_max_pola']  # (A.m-2).
        #           Calculation
        delta_t_pola = delta_t_load_pola + delta_t_break_pola  # s. It is the time of one load.
        tf = delta_t_ini_pola + int(
            i_max_pola / delta_i_pola) * delta_t_pola  # s. It is the polarization current duration.
        n = int(tf / delta_t_pola)  # It is the plot update number.
        time_interval = [0, delta_t_ini_pola + delta_t_pola]  # It is the initial time interval.

        # Dynamic simulation
        for i in range(n):
            Simulator = AlphaPEM(operating_inputs, current_parameters, accessible_physical_parameters,
                                 undetermined_physical_parameters, computing_parameters, initial_variable_values,
                                 time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                t0_interval = Simulator.variables['t'][-1]
                tf_interval = delta_t_ini_pola + (i + 2) * delta_t_pola
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator.solver_variable_names:
                initial_variable_values.append(Simulator.variables[x][-1])

            # Display
            if computing_parameters['type_display'] != "no_display":
                Simulator.Display(ax1, ax2, ax3)

    else:  # elif computing_parameters['type_plot'] == "fixed":
        # Simulation
        Simulator = AlphaPEM(operating_inputs, current_parameters, accessible_physical_parameters,
                             undetermined_physical_parameters, computing_parameters)
        # Display
        if computing_parameters['type_display'] != "no_display":
            Simulator.Display(ax1, ax2, ax3)

    # Plot saving
    Simulator.Save_plot(fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)


def launch_AlphaPEM_for_EIS_current(operating_inputs, current_parameters, accessible_physical_parameters,
                                    undetermined_physical_parameters, computing_parameters):
    """Launch the AlphaPEM simulator for an EIS current density and display the results.

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
            - i0_c_ref : float
                Reference exchange current density at the cathode in A.m-2 (undetermined physical parameter).
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
            - n_gdl : int
                Number of points considered in the GDL (computing parameter).
            - n_mpl : int
                Number of points considered in the MPL (computing parameter).
            - t_purge : tuple
                Time parameters for purging the system (computing parameter).
                It is the purge time interval 'purge_time' and the time between two purges 'delta_purge'.
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
    """

    # Starting time
    start_time = time.time()

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)

    # Initialization
    #       Calculation of the plot update number (n) and the initial time interval (time_interval).
    initial_variable_values = None
    t0_EIS, t_new_start, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = current_parameters['t_EIS']
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = current_parameters['f_EIS']  # These are used for EIS max_step
    #                                                                    actualization.
    f = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)  # It is a list of all the frequency tested.
    n = len(t_new_start)  # It is the plot update number.
    time_interval = [0, t0_EIS]  # It is the initial time interval.

    #       A preliminary simulation run is necessary to equilibrate the internal variables of the cell at i_EIS
    #       prior to initiating the EIS.
    Simulator = AlphaPEM(operating_inputs, current_parameters, accessible_physical_parameters,
                         undetermined_physical_parameters, computing_parameters, initial_variable_values, time_interval)

    # time_interval actualization
    t0_EIS_temp = t0_EIS  # It is the initial time for 1 EIS point.
    tf_EIS_temp = t_new_start[0] + delta_t_break_EIS[0] + delta_t_measurement_EIS[0]  # It is the final time for
    #                                                                                  1 EIS point.
    n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which has been
    #                                                      made.
    time_interval = [t0_EIS_temp, tf_EIS_temp]

    # Recovery of the internal states from the end of the preceding simulation.
    initial_variable_values = []
    for x in Simulator.solver_variable_names:
        initial_variable_values.append(Simulator.variables[x][-1])

    if computing_parameters['type_display'] == "multiple":
        print("A display bug prevents the dynamic updating of the graphs, as it appears that too much data is "
              "involved. However, the data is correctly calculated, and the appropriate plots are saved in the "
              "'results' folder. This display bug does not occur when using a 'synthetic' type_display.")

    # Dynamic simulation
    for i in range(n):
        Simulator = AlphaPEM(operating_inputs, current_parameters, accessible_physical_parameters,
                             undetermined_physical_parameters, computing_parameters, initial_variable_values,
                             time_interval)

        # time_interval actualization
        if i < (n - 1):  # The final simulation does not require actualization.
            t0_EIS_temp = Simulator.variables['t'][-1]  # It is the initial time for 1 EIS point.
            tf_EIS_temp = t_new_start[i + 1] + delta_t_break_EIS[i + 1] + delta_t_measurement_EIS[i + 1]  # It
            #                                                                 is the final time for 1 EIS point.
            n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which
            #                                                      has been made.
            time_interval = [t0_EIS_temp, tf_EIS_temp]  # It is the time interval for 1 EIS point.

        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values = []
        for x in Simulator.solver_variable_names:
            initial_variable_values.append(Simulator.variables[x][-1])

        # Display
        if computing_parameters['type_display'] != "no_display":
            Simulator.Display(ax1, ax2, ax3)

    # Plot saving
    Simulator.Save_plot(fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)
