# ______________________________________Preliminaries_____________________________________

# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib as mpl

# Constants value and functions
from model.AlphaPEM import AlphaPEM
from modules.demo_modules import display_label_operating_inputs_and_physical_parameters, \
    display_value_operating_inputs_and_physical_parameters, \
    recover_for_use_operating_inputs_and_physical_parameters, \
    display_radiobuttons, changeValue, value_control
from configuration.current_densities import step_current, polarization_current, EIS_current

# PyCharm requirement for dynamic plot display
mpl.use("Qt5Agg")


# _______________________________________Functions________________________________________

def create_application():
    """
    This function is responsible for creating the main application window and setting its title. 
    It calls the main_frame() function to create the main graphical elements of the window.
    """
    root = tk.Tk()
    root.title("AlphaPEM")
    # Create the main widget of the application
    main_frame(root)
    root.mainloop()


def main_frame(root):
    """
    This function creates the main graphical elements, such as labels, entry widgets, radio buttons, 
    and buttons, and arranges them in the application window (root). 
    It also initializes the choice dictionary with DoubleVar and IntVar variables for various parameters
    and settings.
    """
    # Creation of the frame
    frame = tk.Frame(root, width=100, height=100)
    frame.grid(row=1, column=0, padx=5, pady=5)

    # Creation of the choice dictionary
    choices = {'Tfc_i': tk.DoubleVar(frame), 'Pa_des_i': tk.DoubleVar(frame), 'Pc_des_i': tk.DoubleVar(frame),
               'Sa_i': tk.DoubleVar(frame), 'Sc_i': tk.DoubleVar(frame), 'Phi_a_des_i': tk.DoubleVar(frame),
               'Phi_c_des_i': tk.DoubleVar(frame), 'Aact_i': tk.DoubleVar(frame), 'Hgdl_i': tk.DoubleVar(frame),
               'Hcl_i': tk.DoubleVar(frame), 'Hmem_i': tk.DoubleVar(frame), 'Hgc_i': tk.DoubleVar(frame),
               'Wgc_i': tk.DoubleVar(frame), 'Lgc_i': tk.DoubleVar(frame), 'epsilon_gdl_i': tk.DoubleVar(frame),
               'epsilon_mc_i': tk.DoubleVar(frame), 'tau_i': tk.DoubleVar(frame), 'epsilon_c_i': tk.DoubleVar(frame),
               'e_i': tk.IntVar(frame), 'Re_i': tk.DoubleVar(frame), 'i0_c_ref_i': tk.DoubleVar(frame),
               'kappa_co_i': tk.DoubleVar(frame), 'kappa_c_i': tk.DoubleVar(frame), 'a_slim_i': tk.DoubleVar(frame),
                'b_slim_i': tk.DoubleVar(frame), 'a_switch_i': tk.DoubleVar(frame), 'C_dl_i': tk.DoubleVar(frame),
               't0_step_i': tk.DoubleVar(frame, 0), 'tf_step_i': tk.DoubleVar(frame, 1000),
               'delta_t_load_step_i': tk.DoubleVar(frame, 20), 'i_ini_step_i': tk.DoubleVar(frame, 0.4),
               'i_final_step_i': tk.DoubleVar(frame, 0.8), 'i_max_pola_i': tk.DoubleVar(frame, 0.5),
               'delta_i_pola_i': tk.DoubleVar(frame, 0.1), 'delta_t_load_pola_i': tk.DoubleVar(frame, 30),
               'delta_t_break_pola_i': tk.DoubleVar(frame, 30),
               'delta_t_ini_pola_i': tk.DoubleVar(frame, 60), 'i_EIS_i': tk.DoubleVar(frame, 0.5),
               'ratio_EIS_i': tk.DoubleVar(frame, 5), 'nb_points_EIS_i': tk.DoubleVar(frame, 50),
               'f_power_min_EIS_i': tk.DoubleVar(frame, -3), 'f_power_max_EIS_i': tk.DoubleVar(frame, 5),
               'nb_f_EIS_i': tk.DoubleVar(frame, 60), 'delta_t_dyn_step_i': tk.DoubleVar(frame, 10),
               't_purge_i': tk.DoubleVar(frame, 0.6), 'delta_t_purge_i': tk.DoubleVar(frame, 15),
               'auxiliaries_choice': tk.IntVar(frame, 0),  'is_purging': tk.IntVar(frame, 0),
               'is_precise': tk.IntVar(frame, 0), 'is_synthetic': tk.IntVar(frame, 0),
               'is_dynamic': tk.IntVar(frame, 0),
               'setting_input': tk.StringVar(frame, 'Enter its specifications')}

    # Displays operating conditions and physical parameters on the screen (without their values)
    display_label_operating_inputs_and_physical_parameters(frame)

    # Displays the value of the operating conditions and physical parameters on the screen.
    label_widgets, entry_widgets = display_value_operating_inputs_and_physical_parameters(frame, choices)

    # Display the radiobuttons on the screen
    display_radiobuttons(frame, choices)

    # Display the 'type of fuel cell' widget on the screen.
    tk.Label(frame, text='Fuel cell:', fg='black', font=('Times New Roman', 12, 'bold')).grid(row=0, column=0,
                                                                                              columnspan=2)
    tk.OptionMenu(frame, choices['setting_input'],
                  'Enter its specifications', 'EH-31 1.5 bar (2021)', 'EH-31 2.0 bar (2021)',
                  'EH-31 2.25 bar (2021)', 'EH-31 2.5 bar (2021)', 'Biao Xie 1.0 bar (2015)',
                  'Biao Xie 1.35 bar (2015)', 'Linhao Fan (2010)',
                  command=lambda value: changeValue(frame, choices, label_widgets, entry_widgets)) \
        .grid(row=0, column=2, columnspan=2)

    # Display the action buttons to select the type of current density to be applied.
    tk.Label(frame, text='Current density:', fg='black', font=('Times New Roman', 12, 'bold')) \
        .grid(row=30, column=0, columnspan=2)
    button_type = {'Step curve': 0, 'Pola curve': 1, 'EIS curve': 2}
    # Button to generate the step curve
    tk.Button(frame, text='Step curve',
              command=lambda: control(choices, button_type['Step curve']), bg='blue', fg='white') \
        .grid(row=30, column=2, padx=10, pady=20)
    # Button to generate the Pola curve
    tk.Button(frame, text='Pola curve', bg='white',
              command=lambda: control(choices, button_type['Pola curve']), fg='black') \
        .grid(row=30, column=3, padx=10, pady=20)
    # Button to generate the EIS curve
    tk.Button(frame, text='EIS curve', bg='brown',
              fg='white', command=lambda: control(choices, button_type['EIS curve'])) \
        .grid(row=30, column=4, padx=10, pady=20)

    # About button
    tk.Button(frame, text='About', bg='black', command=about, fg='white') \
        .grid(row=30, column=5, ipadx=12)


def control(choices, button_type):
    """
    This function is responsible for validating the user inputs by calling the value_control() function. 
    If the input is valid, it then calls the show() function to perform the requested action 
    based on the button_type.
    """
    # Control the values
    value_control(choices)

    # Activate the action
    show(choices, button_type)


def show(choices, button_type):
    """
    This function determines the action to be performed based on the button_type. 
    It calls the AlphaPEM() function with the appropriate parameters and settings to 
    simulate different scenarios such as step curve, polarization curve, or EIS curve.
    """
    if len(choices) == 0:
        return

    # Retrieves parameter values for predefined stacks and keeps them in their standard unit,
    # or converts user-selected quantities into standard units.
    Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, Aact, Hgdl, Hcl, Hmem, Hgc, Wgc, Lgc, epsilon_gdl, epsilon_mc, \
        tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, t_step, i_step, \
        i_max_pola, delta_pola, i_EIS, ratio_EIS, f_EIS, t_EIS, t_purge, delta_t_purge, type_fuel_cell, \
        type_auxiliary, type_purge, max_step, type_display, type_plot \
        = recover_for_use_operating_inputs_and_physical_parameters(choices)

    if button_type == 0:
        type_current = "step_current"
        current_density = step_current  # A.m-2. It is the current density function.
        AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step, i_max_pola,
                 delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau,
                 epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step,
                 n_gdl, t_purge, type_fuel_cell, type_current, type_auxiliary, type_control, type_purge, type_display,
                 type_plot)

    if button_type == 1:
        type_current = "polarization_current"
        current_density = polarization_current  # A.m-2. It is the current density function.
        AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step, i_max_pola,
                 delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau,
                 epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step,
                 n_gdl, t_purge, type_fuel_cell, type_current, type_auxiliary, type_control, type_purge, type_display,
                 type_plot)

    if button_type == 2:
        type_current = "EIS_current"
        current_density = EIS_current  # A.m-2. It is the current density function.
        AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step, i_max_pola,
                 delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau,
                 epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step,
                 n_gdl, t_purge, type_fuel_cell, type_current, type_auxiliary, type_control, type_purge, type_display,
                 type_plot)


def about():
    """
    This function displays information about the program and its author in a dialog box 
    when the "About" button is clicked.
    """
    msg = "AlphaPEM is an open-source PEM fuel cell simulator for control system applications. It is a physical, "\
          "dynamic, two-phase, isothermal, 1D model." \
          "\n\nIt was created by Raphaël GASS, Zhongliang LI, Rachid OUTBIB, Samir JEMEI and Daniel HISSEL." \
          "\n\nIt has been published in the following articles:" \
          "\n    - Gass et al 2024 J. Electrochem. Soc. https://doi.org/10.1149/1945-7111/ad305a," \
          "\n    - Gass et al 2024 SSRN http://dx.doi.org/10.2139/ssrn.4812343." \
          "\n\nContact: gassraphael@proton.me"
    messagebox.showinfo(title='About this program', message=msg)


# ____________________________________Use of the programme ___________________________________
if __name__ == '__main__':
    create_application()
