# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package through a demonstrator interface. Most of the
functionalities are available, but some are not implemented.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import tkinter as tk
from tkinter import messagebox
from ttkthemes import ThemedTk
from tkinter import ttk
import numpy as np
import matplotlib as mpl

# Importing constants' value and functions
from modules.demo_modules import display_label_operating_inputs_and_physical_parameters, \
    display_value_operating_inputs_and_physical_parameters, recover_for_use_operating_inputs_and_physical_parameters, \
    display_radiobuttons, changeValue, value_control, launch_AlphaPEM_for_step_current, \
    launch_AlphaPEM_for_polarization_current, launch_AlphaPEM_for_EIS_current
from configuration.current_densities import step_current, polarization_current, EIS_current

# PyCharm requirement for dynamic plot display
mpl.use("Qt5Agg")


# ____________________________________________________Demo functions____________________________________________________

def create_application():
    """This function creates the main application window and setting its title.
    It calls the main_frame() function to create the main graphical elements of the window.
    """
    root = ThemedTk(theme="arc")
    root.title("AlphaPEM")
    main_frame(root)
    root.mainloop()


def main_frame(root):
    """This function creates the main graphical elements, such as labels, entry widgets, radio buttons, and buttons.
    It arranges them in the application window (root). It also initializes the choice dictionary variables for various
    parameters and settings.

    Parameters:
    -----------
    root : ThemedTk
        The main application window where the graphical elements will be placed.
    """
    # Create the frames
    frame = ttk.Frame(root)
    frame.grid(row=1, column=0, padx=5, pady=5)

    # Create a custom style for the button
    style = ttk.Style()
    style.configure('Blue.TButton', foreground='blue', font=('cmr10', 10, 'bold'))  # Set the font color to blue
    style.configure('Green.TButton', foreground='green', font=('cmr10', 10, 'bold'))  # Set the font color to green
    style.configure('Red.TButton', foreground='red', font=('cmr10', 10, 'bold'))  # Set the font color to red
    style.configure('Black.TButton', foreground='black', font=('cmr10', 10, 'bold'))  # Set the font color to black

    # Create the choice dictionaries
    choices_parameters = {'Tfc (°C)': {'value': tk.DoubleVar(frame), 'label_row': 2, 'label_column': 1},
                          'Pa_des (bar)': {'value': tk.DoubleVar(frame), 'label_row': 2, 'label_column': 3},
                          'Pc_des (bar)': {'value': tk.DoubleVar(frame), 'label_row': 2, 'label_column': 5},
                          'Sa': {'value': tk.DoubleVar(frame), 'label_row': 3, 'label_column': 1},
                          'Sc': {'value': tk.DoubleVar(frame), 'label_row': 3, 'label_column': 3},
                          'Ф_a_des': {'value': tk.DoubleVar(frame), 'label_row': 4, 'label_column': 1},
                          'Ф_c_des': {'value': tk.DoubleVar(frame), 'label_row': 4, 'label_column': 3},
                          'Hgdl (µm)': {'value': tk.DoubleVar(frame), 'label_row': 6, 'label_column': 1},
                          'Hcl (µm)': {'value': tk.DoubleVar(frame), 'label_row': 6, 'label_column': 3},
                          'Hmem (µm)': {'value': tk.DoubleVar(frame), 'label_row': 6, 'label_column': 5},
                          'Hgc (µm)': {'value': tk.DoubleVar(frame), 'label_row': 7, 'label_column': 1},
                          'Wgc (µm)': {'value': tk.DoubleVar(frame), 'label_row': 7, 'label_column': 3},
                          'Lgc (m)': {'value': tk.DoubleVar(frame), 'label_row': 7, 'label_column': 5},
                          'Aact (cm²)': {'value': tk.DoubleVar(frame), 'label_row': 8, 'label_column': 1},
                          'ε_gdl': {'value': tk.DoubleVar(frame, 0.6), 'label_row': 10, 'label_column': 1},
                          'ε_mc': {'value': tk.DoubleVar(frame, 0.25), 'label_row': 10, 'label_column': 3},
                          'τ': {'value': tk.DoubleVar(frame, 1.5), 'label_row': 10, 'label_column': 5},
                          'ε_c': {'value': tk.DoubleVar(frame, 0.2), 'label_row': 11, 'label_column': 1},
                          'e': {'value': tk.IntVar(frame, 4), 'label_row': 11, 'label_column': 3},
                          'Re (µΩ.m²)': {'value': tk.DoubleVar(frame, 1.0), 'label_row': 11, 'label_column': 5},
                          'i0_c_ref (A/m²)': {'value': tk.DoubleVar(frame, 3.0), 'label_row': 12,
                                              'label_column': 1},
                          'κ_co (mol/(m.s.Pa))': {'value': tk.DoubleVar(frame, 1.0), 'label_row': 12,
                                                  'label_column': 3},
                          'κ_c': {'value': tk.DoubleVar(frame, 2.0), 'label_row': 12, 'label_column': 5},
                          'a_slim': {'value': tk.DoubleVar(frame, 0.05), 'label_row': 13, 'label_column': 1},
                          'b_slim': {'value': tk.DoubleVar(frame, 0.1), 'label_row': 13, 'label_column': 3},
                          'a_switch': {'value': tk.DoubleVar(frame, 0.7), 'label_row': 13, 'label_column': 5},
                          'C_dl (MF/m³)': {'value': tk.DoubleVar(frame, 20), 'label_row': 14, 'label_column': 1},
                          't0_step (s)': {'value': tk.DoubleVar(frame, 0), 'label_row': 16, 'label_column': 1},
                          'tf_step (s)': {'value': tk.DoubleVar(frame, 1000), 'label_row': 16, 'label_column': 3},
                          'Δt_load_step (s)': {'value': tk.DoubleVar(frame, 50), 'label_row': 16,
                                               'label_column': 5},
                          'i_ini_step (A/cm²)': {'value': tk.DoubleVar(frame, 0.5), 'label_row': 17,
                                                 'label_column': 1},
                          'i_final_step (A/cm²)': {'value': tk.DoubleVar(frame, 1.5), 'label_row': 17,
                                                   'label_column': 3},
                          'Δt_load_pola (s)': {'value': tk.DoubleVar(frame, 30), 'label_row': 18,
                                               'label_column': 1},
                          'Δt_break_pola (s)': {'value': tk.DoubleVar(frame, 30), 'label_row': 18,
                                                'label_column': 3},
                          'Δt_ini_pola (s)': {'value': tk.DoubleVar(frame, 60), 'label_row': 18,
                                              'label_column': 5},
                          'i_max_pola (A/cm²)': {'value': tk.DoubleVar(frame, 0.5), 'label_row': 19,
                                                 'label_column': 1},
                          'Δi_pola (A/cm²)': {'value': tk.DoubleVar(frame, 0.1), 'label_row': 19,
                                              'label_column': 3},
                          'i_EIS (A/cm²)': {'value': tk.DoubleVar(frame, 0.5), 'label_row': 20,
                                            'label_column': 1},
                          'ratio_EIS (%)': {'value': tk.DoubleVar(frame, 5), 'label_row': 20, 'label_column': 3},
                          'nb_points_EIS': {'value': tk.IntVar(frame, 50), 'label_row': 20, 'label_column': 5},
                          'f_power_min_EIS': {'value': tk.IntVar(frame, -3), 'label_row': 21, 'label_column': 1},
                          'f_power_max_EIS': {'value': tk.IntVar(frame, 5), 'label_row': 21, 'label_column': 3},
                          'nb_f_EIS': {'value': tk.IntVar(frame, 60), 'label_row': 21, 'label_column': 5},
                          'Δt_dyn_step (s)': {'value': tk.DoubleVar(frame, 10), 'label_row': 23,
                                              'label_column': 1},
                          't_purge (s)': {'value': tk.DoubleVar(frame, 0.6), 'label_row': 23, 'label_column': 3},
                          'Δt_purge (s)': {'value': tk.DoubleVar(frame, 15), 'label_row': 23, 'label_column': 5},
                          'max_step (s)': {'value': tk.DoubleVar(frame, 0.1), 'label_row': 24, 'label_column': 1},
                          'n_gdl': {'value': tk.IntVar(frame, 10), 'label_row': 24, 'label_column': 3}}

    choices_buttons = {'type_fuel_cell': {'value': tk.StringVar(frame, 'Enter your specifications'),
                                          'label_row': 0},
                       'type_auxiliary': {'value': tk.IntVar(frame, 2), 'label_row': 26},
                       'type_control': {'value': tk.IntVar(frame, 0), 'label_row': 27},
                       'type_purge': {'value': tk.IntVar(frame, 0), 'label_row': 28},
                       'type_display': {'value': tk.IntVar(frame, 1), 'label_row': 29},
                       'type_plot': {'value': tk.IntVar(frame, 0), 'label_row': 30}}

    # Displays operating conditions and physical parameters on the screen (without their values)
    display_label_operating_inputs_and_physical_parameters(frame, choices_parameters)

    # Displays the value of the operating conditions and physical parameters on the screen.
    label_widgets, entry_widgets = display_value_operating_inputs_and_physical_parameters(frame, choices_parameters)

    # Display the radiobuttons on the screen
    display_radiobuttons(frame, choices_buttons)

    # Display the 'type of fuel cell' widget on the screen.
    ttk.Label(frame, text='Fuel cell:', font=('cmr10', 12, 'bold')).grid(row=0, column=0, columnspan=2)
    ttk.OptionMenu(frame, choices_buttons['type_fuel_cell']['value'],
                   'Enter your specifications', 'Enter your specifications', 'EH-31 1.5 bar (2021)',
                   'EH-31 2.0 bar (2021)', 'EH-31 2.25 bar (2021)', 'EH-31 2.5 bar (2021)', 'Biao Xie 1.0 bar (2015)',
                   'Biao Xie 1.35 bar (2015)', 'Linhao Fan (2010)',
                   command=lambda value: changeValue(frame, choices_parameters, choices_buttons, label_widgets,
                                                     entry_widgets)) \
                                                                   .grid(row=0, column=2, columnspan=2)

    # Display the action buttons to select the type of current density to be applied.
    ttk.Label(frame, text='Current density:', font=('cmr10', 12, 'bold')).grid(row=31, column=0, columnspan=2)
    current_button = {'Step curve': 0, 'Pola curve': 1, 'EIS curve': 2}
    #       Button to generate the step curve
    ttk.Button(frame, text='Step curve', style='Blue.TButton',
               command=lambda: control_current_button(choices_parameters, choices_buttons,
                                                      current_button['Step curve'])) \
                                                                         .grid(row=31, column=2, padx=10, pady=20)
    #       Button to generate the Pola curve
    ttk.Button(frame, text='Pola curve', style='Green.TButton',
               command=lambda: control_current_button(choices_parameters, choices_buttons,
                                                      current_button['Pola curve'])) \
                                                                         .grid(row=31, column=3, padx=10, pady=20)
    #       Button to generate the EIS curve
    ttk.Button(frame, text='EIS curve', style='Red.TButton',
               command=lambda: control_current_button(choices_parameters, choices_buttons,
                                                      current_button['EIS curve'])) \
                                                                         .grid(row=31, column=4, padx=10, pady=20)
    #       About button
    ttk.Button(frame, text='About', style='Black.TButton', command=about) \
                                                                         .grid(row=31, column=5, ipadx=12)


def control_current_button(choices_parameters, choices_buttons, current_button):
    """This function is responsible for validating the user inputs by calling the value_control() function. If the input
    is valid, it then calls the show_current_button function to perform the requested action based on the button_type.

    Parameters
    ----------
    choices_parameters : dict
        A dictionary containing the parameter information.
    choices_buttons : dict
        A dictionary containing the button information.
    current_button : dict
        A dictionary representing the clicked button.
    """
    # Control the values
    value_control(choices_parameters, choices_buttons, current_button)

    # Activate the action
    show_current_button(choices_parameters, choices_buttons, current_button)


def show_current_button(choices_parameters, choices_buttons, current_button):
    """This function determines the action to be performed based on the button_type.

    Parameters
    ----------
    choices_parameters : dict
        A dictionary containing the parameter information.
    choices_buttons : dict
        A dictionary containing the button information.
    current_button : dict
        A dictionary representing the clicked button.
    """

    # Retrieves parameter values for predefined stacks and keeps them in their standard unit, or converts user-selected
    # quantities into standard units.
    Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, Aact, Hgdl, Hcl, Hmem, Hgc, Wgc, Lgc, epsilon_gdl, \
        epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, t_step, \
        i_step, i_max_pola, delta_pola, i_EIS, ratio_EIS, f_EIS, t_EIS, t_purge, delta_t_purge, max_step, n_gdl, \
        type_fuel_cell, type_auxiliary, type_control, type_purge, type_display, type_plot \
        = recover_for_use_operating_inputs_and_physical_parameters(choices_parameters, choices_buttons)

    if current_button == 0:
        type_current = "step"
        current_density = step_current
        launch_AlphaPEM_for_step_current(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step,
                                         i_step, i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                         Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re,
                                         i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl,
                                         t_purge, type_fuel_cell, type_current, type_auxiliary, type_control,
                                         type_purge, type_display, type_plot)

    if current_button == 1:
        type_current = "polarization"
        current_density = polarization_current
        launch_AlphaPEM_for_polarization_current(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des,
                                                 t_step, i_step, i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS,
                                                 Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc,
                                                 epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                                 a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell, type_current,
                                                 type_auxiliary, type_control, type_purge, type_display, type_plot)

    if current_button == 2:
        type_current = "EIS"
        current_density = EIS_current
        launch_AlphaPEM_for_EIS_current(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step,
                                        i_step, i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                        Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re,
                                        i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl,
                                        t_purge, type_fuel_cell, type_current, type_auxiliary, type_control, type_purge,
                                        type_display, type_plot)


def about():
    """This function displays information about the program and its author in a dialog box when the "About" button is
    clicked.
    """
    msg = "AlphaPEM is an open-source PEM fuel cell simulator for control system applications. It is a physical, " \
          "dynamic, two-phase, isothermal, 1D model." \
          "\n\nIt was created by Raphaël GASS, Zhongliang LI, Rachid OUTBIB, Samir JEMEI and Daniel HISSEL." \
          "\n\nIt has been published in the following articles:" \
          "\n    - Gass et al 2024 J. Electrochem. Soc. https://doi.org/10.1149/1945-7111/ad305a," \
          "\n    - Gass et al 2024 SSRN http://dx.doi.org/10.2139/ssrn.4812343." \
          "\n\nContact: gassraphael@proton.me"
    messagebox.showinfo(title='About this program', message=msg)


# ________________________________________________Use of the programme _________________________________________________
if __name__ == '__main__':
    create_application()
