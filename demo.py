# ______________________________________Preliminaries_____________________________________

# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import messagebox
import numpy as np

# Constants value and functions
from AlphaPEM_main import AlphaPEM
from modules.interface_modules import display_label_operating_inputs_and_physical_parameters, \
    display_value_operating_inputs_and_physical_parameters, \
    recover_for_use_operating_inputs_and_physical_parameters, \
    display_radiobuttons, changeValue, value_control
from configuration.current_densities import step_current, polarization_current, EIS_current
from modules.settings_modules import EIS_parameters

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

    # Creation of the choice dictionnary
    choices = {'Tfc_i': tk.DoubleVar(frame),
               'Pa_des_i': tk.DoubleVar(frame), 'Pc_des_i': tk.DoubleVar(frame),
               'Sa_i': tk.DoubleVar(frame), 'Sc_i': tk.DoubleVar(frame),
               'Phi_a_des_i': tk.DoubleVar(frame), 'Phi_c_des_i': tk.DoubleVar(frame),
               'Aact_i': tk.DoubleVar(frame),
               'Hgdl_i': tk.DoubleVar(frame), 'Hcl_i': tk.DoubleVar(frame), 'Hmem_i': tk.DoubleVar(frame),
               'Hgc_i': tk.DoubleVar(frame), 'Wgc_i': tk.DoubleVar(frame), 'Lgc_i': tk.DoubleVar(frame),
               'epsilon_gdl_i': tk.DoubleVar(frame), 'epsilon_cl_i': tk.DoubleVar(frame),
               'epsilon_mc_i': tk.DoubleVar(frame),
               'tau_i': tk.DoubleVar(frame), 'epsilon_c_i': tk.DoubleVar(frame), 'e_i': tk.IntVar(frame),
               'Re_i': tk.DoubleVar(frame), 'i0_c_ref_i': tk.DoubleVar(frame), 'kappa_co_i': tk.DoubleVar(frame),
               'kappa_c_i': tk.DoubleVar(frame),
               't0_i': tk.DoubleVar(frame, 0), 'tf_i': tk.DoubleVar(frame, 500),
               'delta_t_load_i': tk.DoubleVar(frame, 20),
               'i_ini_i': tk.DoubleVar(frame, 0.4), 'i_final_i': tk.DoubleVar(frame, 0.8),
               'i_pola_i': tk.DoubleVar(frame, 0.5), 'i_EIS_i': tk.DoubleVar(frame, 0.5),
               't_purge_i': tk.DoubleVar(frame, 0.6), 'delta_t_purge_i': tk.DoubleVar(frame, 15),
               'auxiliaries_choice': tk.IntVar(frame, 0), 'is_purging': tk.IntVar(frame, 0),
               'is_precise': tk.IntVar(frame, 0), 'is_synthetic': tk.IntVar(frame, 0),
               'is_dynamic': tk.IntVar(frame, 0), 'setting_input': tk.StringVar(frame, 'Write them')}

    # Displays operating conditions and physical parameters on the screen (without their values)
    display_label_operating_inputs_and_physical_parameters(frame)

    # Displays the value of the operating conditions and physical parameters on the screen.
    label_widgets, entry_widgets = display_value_operating_inputs_and_physical_parameters(frame, choices)

    # Display the radiobuttons on the screen
    display_radiobuttons(frame, choices)

    # Display the setting widget on the screen.
    tk.Label(frame, text='Settings choice:', fg='black', font=('Times New Roman', 12)).grid(row=0, column=1)
    tk.OptionMenu(frame, choices['setting_input'],
                  'Write them', 'Raphaël Gass 1', 'Biao Xie 1(2015)', 'Linhao Fan',
                  command=lambda value: changeValue(frame, choices, label_widgets, entry_widgets)) \
        .grid(row=0, column=2, columnspan=2)

    # Display the action buttons to click to launch the program according to what you want to do
    button_type = {'Current response': 0, 'Polarisation curve': 1, 'EIS curve (WIP)': 2}
    # Button to generate the current response curve
    tk.Button(frame, text='Current response',
              command=lambda: control(choices, button_type['Current response']), bg='blue', fg='white') \
        .grid(row=25, column=0, padx=10, pady=20)
    # Button to generate the polarisation curve
    tk.Button(frame, text='Polarisation curve', bg='white',
              command=lambda: control(choices, button_type['Polarisation curve']), fg='black') \
        .grid(row=25, column=1, padx=10, pady=20)
    # Button to generate the EIS curve
    tk.Button(frame, text='EIS curve (WIP)', bg='brown',
              fg='white', command=lambda: control(choices, button_type['EIS curve (WIP)'])) \
        .grid(row=25, column=2, padx=10, pady=20)

    # About button
    tk.Button(frame, text='About', bg='black', command=about, fg='white') \
        .grid(row=25, column=3, ipadx=12)


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
    simulate different scenarios such as current response, polarization curve, or EIS curve.
    """
    if len(choices) == 0:
        return

    # Retrieves parameter values for predefined stacks and keeps them in their standard unit,
    # or converts user-selected quantities into standard units.
    Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, \
        Aact, Hgdl, Hcl, Hmem, Hgc, Wgc, Lgc, \
        epsilon_gdl, epsilon_cl, epsilon_mc, tau, epsilon_c, e, \
        Re, i0_c_ref, kappa_co, kappa_c, \
        t_step, i_step, i_pola, i_EIS, t_purge, delta_t_purge, \
        type_fuel_cell, type_auxiliary, type_purge, max_step, type_display, type_plot \
        = recover_for_use_operating_inputs_and_physical_parameters(choices)

    if button_type == 0:
        type_current = "step_current"
        delta_pola = np.nan, np.nan, np.nan  # It is the loading time, breaking time,
        # and current density variation at each time step.
        t_EIS, f_EIS = np.nan, np.nan  # It is the EIS parameters.
        delta_t_plot = 1  # s. It is the interval between two dynamic plots update.
        current_density = step_current  # A.m-2. It is the current density function.
        AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des,
                 t_step, i_step, i_pola, delta_pola, i_EIS, t_EIS, f_EIS, t_purge,
                 Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc,
                 epsilon_gdl, epsilon_cl, tau, epsilon_mc, epsilon_c, e, kappa_co, Re, i0_c_ref, kappa_c,
                 delta_t_plot, max_step, type_fuel_cell, type_current, type_auxiliary, type_purge, type_display, type_plot)

    if button_type == 1:
        type_current = "polarization_current"
        delta_pola = 30, 5 * 60, 100  # s, s, A.m-2. It is the loading time, breaking time,
        # and current density variation at each time step.
        t_EIS, f_EIS = np.nan, np.nan  # It is the EIS parameters.
        delta_t_plot = 1  # s. It is the interval between two dynamic plots update.
        current_density = polarization_current  # A.m-2. It is the current density function.
        AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des,
                 t_step, i_step, i_pola, delta_pola, i_EIS, t_EIS, f_EIS, t_purge,
                 Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc,
                 epsilon_gdl, epsilon_cl, tau, epsilon_mc, epsilon_c, e, kappa_co, Re, i0_c_ref, kappa_c,
                 delta_t_plot, max_step, type_fuel_cell, type_current, type_auxiliary, type_purge, type_display, type_plot)

    if button_type == 2:
        type_current = "EIS_current"
        delta_pola = np.nan, np.nan, np.nan  # It is the loading time, breaking time,
        # and current density variation at each time step.
        t_EIS, f_EIS = EIS_parameters()  # (s, s, s), (s-1, s-1, ). It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        delta_t_plot = 1  # s. It is the interval between two dynamic plots update.
        current_density = EIS_current  # A.m-2. It is the current density function.
        AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des,
                 t_step, i_step, i_pola, delta_pola, i_EIS, t_EIS, f_EIS, t_purge,
                 Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc,
                 epsilon_gdl, epsilon_cl, tau, epsilon_mc, epsilon_c, e, kappa_co, Re, i0_c_ref, kappa_c,
                 delta_t_plot, max_step, type_fuel_cell, type_current, type_auxiliary, type_purge, type_display, type_plot)


def about():
    """
    This function displays information about the program and its author in a dialog box 
    when the "About" button is clicked.
    """
    msg = "AlphaPEM models the PEM fuel cell in a 9-nodes, dynamic, two-phase, 1D model." \
          "\nIt was created by Raphaël Gass." \
          "\nContact: gassraphael@gmail.com"
    messagebox.showinfo(title='About this program', message=msg)


# ____________________________________Use of the programme ___________________________________
if __name__ == '__main__':
    create_application()
