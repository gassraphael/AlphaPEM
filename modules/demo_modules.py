# -*- coding: utf-8 -*-
import numpy as np
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk

# Constants value and functions
from model.AlphaPEM import AlphaPEM
from modules.settings_modules import stored_operating_inputs, stored_physical_parameters, EIS_parameters
from modules.main_modules import figures_preparation, plot_saving

"""
This modul contains some of the required functions for the main program AlphaPEM_interface.
"""


def changeValue(frame, choices_parameters, choices_buttons, Label_widgets, Entry_widgets):
    """This function is called when the user selects a specific option from a dropdown menu
    (choices_buttons['type_fuel_cell']). Depending on the selected option, it either hides or displays specific
    input fields (labels or entry widgets) on the GUI.
    """

    if choices_buttons['type_fuel_cell']['value'].get() != 'Enter your specifications':
        # Hide all the Entries if they are already displayed
        for k, v in Entry_widgets.items():
            v.grid_forget()

        # Recovers the new settings
        recover_for_display_operating_inputs_and_physical_parameters(frame, choices_parameters, choices_buttons)

        # Display the Labels
        for k, v in choices_parameters.items():
            ttk.Label(frame, width=7, anchor='w', textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)

    else:  # choices_buttons['type_fuel_cell']['value'].get() == 'Enter your specifications':
        # Hide the Labels if they are already displayed
        for k, v in Label_widgets.items():
            v.grid_forget()

        # Saves and displays the user entries
        for k, v in choices_parameters.items():
            ttk.Label(frame, width=7, anchor='w', textvariable=v['value']). \
                grid(row=v['label_row'], column=v['label_column'], padx=5)

def display_label_operating_inputs_and_physical_parameters(frame, choices_parameters):
    """
    This function displays labels on the GUI, representing operating conditions and physical parameters, 
    without their actual values.
    """

    # Display the titles
    ttk.Label(frame, text='Operating conditions', font=('cmr10', 12, 'bold')). \
        grid(row=1, column=0, columnspan=6, ipady=15)
    ttk.Label(frame, text='Accessible physical parameters', font=('cmr10', 12, 'bold')). \
        grid(row=5, column=0, columnspan=6, ipady=15)
    ttk.Label(frame, text='Undetermined physical parameters', font=('cmr10', 12, 'bold')). \
        grid(row=9, column=0, columnspan=6, ipady=15)
    ttk.Label(frame, text='Current density parameters', font=('cmr10', 12, 'bold')). \
        grid(row=15, column=0, columnspan=6, ipady=15)
    ttk.Label(frame, text='Computing parameters', font=('cmr10', 12, 'bold')). \
        grid(row=22, column=0, columnspan=6, ipady=15)

    # Display the parameters labels
    for k, p in choices_parameters.items():
        ttk.Label(frame, text=k, font=('cmr10', 10)).grid(row=p['label_row'], column=p['label_column'] - 1, sticky="w")


def display_value_operating_inputs_and_physical_parameters(frame, choices_parameters):
    """This function displays entry widgets on the GUI,
    where the user can enter values for operating conditions and physical parameters.
    """
    Label_widgets = {}
    for k, p in choices_parameters.items():
        Label_widgets['Label ' + k] = ttk.Label(frame, width=7, anchor='w', textvariable=p['value'])

    Entry_widgets = {}
    for k, p in choices_parameters.items():
        Entry_widgets['Entry ' + k] = ttk.Entry(frame, width=7, textvariable=p['value'])
        Entry_widgets['Entry ' + k].grid(row=p['label_row'], column=p['label_column'], padx=5)

    return Label_widgets, Entry_widgets


def display_radiobuttons(frame, choices_buttons):
    """
    This function displays radiobuttons on the GUI, allowing the user to make choices for purging, 
    results display, plot style, etc.
    """
    ttk.Label(frame, text='Model possibilities', font=('cmr10', 12, 'bold')) \
        .grid(row=25, column=0, columnspan=6, ipady=15)

    # Ask the user to choose an option and save it
    ttk.Label(frame, text='Auxiliaries: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=0, sticky="w")
    ttk.Radiobutton(frame, text='No auxiliaries', value=0, variable=choices_buttons['type_auxiliary']['value']). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=1, sticky="w")
    ttk.Radiobutton(frame, text='Closed anode', value=1, variable=choices_buttons['type_auxiliary']['value']). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(frame, text='Opened anode', value=2, variable=choices_buttons['type_auxiliary']['value']). \
        grid(row=choices_buttons['type_auxiliary']['label_row'], column=3, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(frame, text='Control: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_control']['label_row'], column=0, sticky="w")
    ttk.Radiobutton(frame, text='No control', value=0, variable=choices_buttons['type_control']['value']). \
        grid(row=choices_buttons['type_control']['label_row'], column=1, sticky="w")
    ttk.Radiobutton(frame, text='Humidity', value=1, variable=choices_buttons['type_control']['value']). \
        grid(row=choices_buttons['type_control']['label_row'], column=2, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(frame, text='Purge: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_purge']['label_row'], column=0, sticky="w")
    ttk.Radiobutton(frame, text='No purge', value=0, variable=choices_buttons['type_purge']['value']). \
        grid(row=choices_buttons['type_purge']['label_row'], column=1, sticky="w")
    ttk.Radiobutton(frame, text='Periodic', value=1, variable=choices_buttons['type_purge']['value']). \
        grid(row=choices_buttons['type_purge']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(frame, text='Constant', value=2, variable=choices_buttons['type_purge']['value']). \
        grid(row=choices_buttons['type_purge']['label_row'], column=3, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(frame, text='Display: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_display']['label_row'], column=0, sticky="w")
    ttk.Radiobutton(frame, text='No display', value=0, variable=choices_buttons['type_display']['value']). \
        grid(row=choices_buttons['type_display']['label_row'], column=1, sticky="w")
    ttk.Radiobutton(frame, text='Synthetic', value=1, variable=choices_buttons['type_display']['value']). \
        grid(row=choices_buttons['type_display']['label_row'], column=2, sticky="w")
    ttk.Radiobutton(frame, text='Multiple', value=2, variable=choices_buttons['type_display']['value']). \
        grid(row=choices_buttons['type_display']['label_row'], column=3, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(frame, text='Plot: ', font=('cmr10', 12)). \
        grid(row=choices_buttons['type_plot']['label_row'], column=0, sticky="w")
    ttk.Radiobutton(frame, text='Fixed', value=0, variable=choices_buttons['type_plot']['value']). \
        grid(row=choices_buttons['type_plot']['label_row'], column=1, sticky="w")
    ttk.Radiobutton(frame, text='Dynamic', value=1, variable=choices_buttons['type_plot']['value']). \
        grid(row=choices_buttons['type_plot']['label_row'], column=2, sticky="w")


def recover_for_display_operating_inputs_and_physical_parameters(frame, choices_parameters, choices_buttons):
    """
    This function retrieves parameter values for predefined stacks (e.g., "EH-31 1.5 bar (2021)", "Biao Xie 1.0 bar
    (2015)", etc.) and converts them to appropriate units for display on the GUI.
    """

    if choices_buttons['type_fuel_cell']['value'].get() == "EH-31 1.5 bar (2021)": type_fuel_cell = "EH-31_1.5"
    elif choices_buttons['type_fuel_cell']['value'].get() == "EH-31 2.0 bar (2021)": type_fuel_cell = "EH-31_2.0"
    elif choices_buttons['type_fuel_cell']['value'].get() == "EH-31 2.25 bar (2021)": type_fuel_cell = "EH-31_2.25"
    elif choices_buttons['type_fuel_cell']['value'].get() == "EH-31 2.5 bar (2021)": type_fuel_cell = "EH-31_2.5"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Biao Xie 1.0 bar (2015)": type_fuel_cell = "BX_1.0"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Biao Xie 1.35 bar (2015)": type_fuel_cell = "BX_1.35"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Linhao Fan (2010)": type_fuel_cell = "LF"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Enter your specifications": type_fuel_cell = "manual_setup"
    else: raise ValueError('the type_fuel_cell given is not valid.')

    Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, i_max_pola = stored_operating_inputs(type_fuel_cell)

    Hcl, epsilon_mc, tau, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, Re, i0_c_ref, kappa_co, kappa_c, \
        a_slim, b_slim, a_switch, C_dl = \
        stored_physical_parameters(type_fuel_cell)

    choices_parameters['Tfc (°C)']['value'].set(np.round(Tfc - 273.15))  # °C
    choices_parameters['Pa_des (bar)']['value'].set(np.round(Pa_des / 1e5, 2))  # bar
    choices_parameters['Pc_des (bar)']['value'].set(np.round(Pc_des / 1e5, 2))  # bar
    choices_parameters['Sa']['value'].set(np.round(Sa, 1))
    choices_parameters['Sc']['value'].set(np.round(Sc, 1))
    choices_parameters['Ф_a_des']['value'].set(np.round(Phi_a_des, 1))
    choices_parameters['Ф_c_des']['value'].set(np.round(Phi_c_des, 1))
    choices_parameters['Aact (cm²)']['value'].set(np.round(Aact * 1e4))  # cm²
    choices_parameters['Hgdl (µm)']['value'].set(np.round(Hgdl * 1e6))  # µm
    choices_parameters['Hcl (µm)']['value'].set(np.round(Hcl * 1e6))  # µm
    choices_parameters['Hmem (µm)']['value'].set(np.round(Hmem * 1e6))  # µm
    choices_parameters['Hgc (µm)']['value'].set(np.round(Hgc * 1e6))  # µm
    choices_parameters['Wgc (µm)']['value'].set(np.round(Wgc * 1e6))  # µm
    choices_parameters['Lgc (m)']['value'].set(np.round(Lgc, 2))  # µm
    choices_parameters['ε_gdl']['value'].set(np.round(epsilon_gdl, 2))
    choices_parameters['ε_mc']['value'].set(np.round(epsilon_mc, 2))
    choices_parameters['τ']['value'].set(np.round(tau, 2))
    choices_parameters['ε_c']['value'].set(np.round(epsilon_c, 2))
    choices_parameters['e']['value'].set(e)
    choices_parameters['Re (µΩ.m²)']['value'].set(np.round(Re * 1e6, 2))  # µΩ.m²
    choices_parameters['i0_c_ref (A/m²)']['value'].set(np.round(i0_c_ref, 2))  # A.m-2
    choices_parameters['κ_co (mol/(m.s.Pa))']['value'].set(np.round(kappa_co, 2))  # mol.m-1.s-1.Pa-1
    choices_parameters['κ_c']['value'].set(np.round(kappa_c, 2))
    choices_parameters['a_slim']['value'].set(np.round(a_slim, 4))
    choices_parameters['b_slim']['value'].set(np.round(b_slim, 4))
    choices_parameters['a_switch']['value'].set(np.round(a_switch, 4))
    choices_parameters['C_dl (MF/m³)']['value'].set(np.round(C_dl * 1e-6, 2))  # MF.m-3
    choices_parameters['i_max_pola (A/cm²)']['value'].set(np.round(i_max_pola / 1e4, 2))  # A/cm²


def recover_for_use_operating_inputs_and_physical_parameters(choices_parameters, choices_buttons):
    """This function retrieves and converts the parameter values from the GUI into standard units
    for further calculations.
    """
    Tfc = choices_parameters['Tfc (°C)']['value'].get() + 273.15  # K
    Pa_des = choices_parameters['Pa_des (bar)']['value'].get() * 1e5  # Pa
    Pc_des = choices_parameters['Pc_des (bar)']['value'].get() * 1e5  # Pa
    Sa = choices_parameters['Sa']['value'].get()
    Sc = choices_parameters['Sc']['value'].get()
    Phi_a_des = choices_parameters['Ф_a_des']['value'].get()
    Phi_c_des = choices_parameters['Ф_c_des']['value'].get()
    Aact = choices_parameters['Aact (cm²)']['value'].get() * 1e-4  # m²
    Hgdl = choices_parameters['Hgdl (µm)']['value'].get() * 1e-6  # m
    Hcl = choices_parameters['Hcl (µm)']['value'].get() * 1e-6  # m
    Hmem = choices_parameters['Hmem (µm)']['value'].get() * 1e-6  # m
    Hgc = choices_parameters['Hgc (µm)']['value'].get() * 1e-6  # m
    Wgc = choices_parameters['Wgc (µm)']['value'].get() * 1e-6  # m
    Lgc = choices_parameters['Lgc (m)']['value'].get()  # m
    epsilon_gdl = choices_parameters['ε_gdl']['value'].get()
    epsilon_mc = choices_parameters['ε_mc']['value'].get()
    tau = choices_parameters['τ']['value'].get()
    epsilon_c = choices_parameters['ε_c']['value'].get()
    e = choices_parameters['e']['value'].get()
    Re = choices_parameters['Re (µΩ.m²)']['value'].get() * 1e-6  # ohm.m²
    i0_c_ref = choices_parameters['i0_c_ref (A/m²)']['value'].get()  # A.m-2
    kappa_co = choices_parameters['κ_co (mol/(m.s.Pa))']['value'].get()  # mol.m-1.s-1.Pa-1
    kappa_c = choices_parameters['κ_c']['value'].get()
    a_slim = choices_parameters['a_slim']['value'].get()
    b_slim = choices_parameters['b_slim']['value'].get()
    a_switch = choices_parameters['a_switch']['value'].get()
    C_dl = choices_parameters['C_dl (MF/m³)']['value'].get() * 1e6  # F.m-3
    t_step = (choices_parameters['t0_step (s)']['value'].get(), choices_parameters['tf_step (s)']['value'].get(),
              choices_parameters['Δt_load_step (s)']['value'].get(),
              choices_parameters['Δt_dyn_step (s)']['value'].get())  # (s, s, s, s)
    i_step = (choices_parameters['i_ini_step (A/cm²)']['value'].get() * 1e4,
              choices_parameters['i_final_step (A/cm²)']['value'].get() * 1e4)  # (A.m-2, A.m-2)
    i_max_pola = choices_parameters['i_max_pola (A/cm²)']['value'].get() * 1e4  # A.m-2
    delta_pola = (choices_parameters['Δt_load_pola (s)']['value'].get(),
                  choices_parameters['Δt_break_pola (s)']['value'].get(),
                  choices_parameters['Δi_pola (A/cm²)']['value'].get() * 1e4,
                  choices_parameters['Δt_ini_pola (s)']['value'].get()) # (s, s, A.m-2, s)
    i_EIS = choices_parameters['i_EIS (A/cm²)']['value'].get() * 1e4  # (A.m-2)
    ratio_EIS = choices_parameters['ratio_EIS (%)']['value'].get() / 100
    f_EIS = (choices_parameters['f_power_min_EIS']['value'].get(),
             choices_parameters['f_power_max_EIS']['value'].get(),
             choices_parameters['nb_f_EIS']['value'].get(),
             choices_parameters['nb_points_EIS']['value'].get())
    t_EIS = EIS_parameters(f_EIS)  # Time parameters for the EIS_current density function.
    t_purge = choices_parameters['t_purge (s)']['value'].get()  # s
    delta_t_purge = choices_parameters['Δt_purge (s)']['value'].get()  # s
    max_step = choices_parameters['max_step (s)']['value'].get()  # s
    n_gdl = choices_parameters['n_gdl']['value'].get()

    if choices_buttons['type_fuel_cell']['value'].get() == "EH-31 1.5 bar (2021)":
        type_fuel_cell = "EH-31_1.5"
    elif choices_buttons['type_fuel_cell']['value'].get() == "EH-31 2.0 bar (2021)":
        type_fuel_cell = "EH-31_2.0"
    elif choices_buttons['type_fuel_cell']['value'].get() == "EH-31 2.25 bar (2021)":
        type_fuel_cell = "EH-31_2.25"
    elif choices_buttons['type_fuel_cell']['value'].get() == "EH-31 2.5 bar (2021)":
        type_fuel_cell = "EH-31_2.5"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Biao Xie 1.0 bar (2015)":
        type_fuel_cell = "BX_1.0"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Biao Xie 1.35 bar (2015)":
        type_fuel_cell = "BX_1.35"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Linhao Fan (2010)":
        type_fuel_cell = "LF"
    elif choices_buttons['type_fuel_cell']['value'].get() == "Enter your specifications":
        type_fuel_cell = "manual_setup"

    if choices_buttons['type_auxiliary']['value'].get() == 0:
        type_auxiliary = "no_auxiliary"
    elif choices_buttons['type_auxiliary']['value'].get() == 1:
        type_auxiliary = "closed_anode"
    else:
        type_auxiliary = "opened_anode"

    if choices_buttons['type_control']['value'].get() == 0:
        type_control = "no_control"
    else:
        type_control = "Phi_des"

    if choices_buttons['type_purge']['value'].get() == 0:
        type_purge = "no_purge"
    elif choices_buttons['type_purge']['value'].get() == 1:
        type_purge = "periodic_purge"
    else:
        type_purge = "constant_purge"

    if choices_buttons['type_display']['value'].get() == 0:
        type_display = "no_display"
    elif choices_buttons['type_display']['value'].get() == 1:
        type_display = "synthetic"
    else:
        type_display = "multiple_display"

    if choices_buttons['type_plot']['value'].get() == 0:
        type_plot = "fixed"
    else:
        type_plot = "dynamic"

    return (Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, Aact, Hgdl, Hcl, Hmem, Hgc, Wgc, Lgc, epsilon_gdl,
            epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, t_step,
            i_step, i_max_pola, delta_pola, i_EIS, ratio_EIS, f_EIS, t_EIS, t_purge, delta_t_purge, max_step, n_gdl,
            type_fuel_cell, type_auxiliary, type_control, type_purge, type_display, type_plot)


def value_control(choices_parameters, choices_buttons, current_button):
    """
    This function checks the integrity of the values entered by the user and 
    returns an empty tuple if they are not valid.
    """

    # The values entered by the user are checked for compliance
    if choices_parameters['Tfc (°C)']['value'].get() < 0:
        messagebox.showerror(title='Temperatures', message='Negative temperatures do not exist in the Kelvin scale.')
        choices.clear()
        return
    if choices_parameters['Pa_des (bar)']['value'].get() < 0 or \
            choices_parameters['Pc_des (bar)']['value'].get() < 0 or \
            choices_parameters['Pc_des (bar)']['value'].get() > 5.0 or \
            choices_parameters['Pc_des (bar)']['value'].get() > 5.0:
        messagebox.showerror(title='Desired pressures', message=
        'Desired pressure should be positive and bellow 5.0 bars.')
        choices.clear()
        return
    if choices_parameters['Sa']['value'].get() < 1 or choices_parameters['Sa']['value'].get() > 5 or \
            choices_parameters['Sc']['value'].get() < 1 or choices_parameters['Sc']['value'].get() > 5:
        messagebox.showerror(title='Stoichiometric ratios', message=
        'The stoichiometric ratios Sa and Sc should be between 1 and 5.')
        choices.clear()
        return
    if choices_parameters['Ф_a_des']['value'].get() < 0 or choices_parameters['Ф_a_des']['value'].get() > 1 or \
            choices_parameters['Ф_c_des']['value'].get() < 0 or \
            choices_parameters['Ф_c_des']['value'].get() > 1:
        messagebox.showerror(title='Desired humidity', message='The desired humidities should be between 0 and 1.')
        choices.clear()
        return
    if choices_parameters['Aact (cm²)']['value'].get() < 0:
        messagebox.showerror(title='Active area', message='Negative active area is impossible.')
        choices.clear()
        return
    if choices_parameters['Hgdl (µm)']['value'].get() < 1 or \
            choices_parameters['Hgdl (µm)']['value'].get() > 1000 or \
            choices_parameters['Hcl (µm)']['value'].get() < 1 or \
            choices_parameters['Hcl (µm)']['value'].get() > 1000 or \
            choices_parameters['Hmem (µm)']['value'].get() < 1 or \
            choices_parameters['Hmem (µm)']['value'].get() > 1000:
        messagebox.showerror(title='MEA thickness', message=
        'All MEA components generally have a thickness between 1µm and 1mm.')
        choices.clear()
        return
    if choices_parameters['Hgc (µm)']['value'].get() < 10 or \
            choices_parameters['Hgc (µm)']['value'].get() > 10000 or \
            choices_parameters['Wgc (µm)']['value'].get() < 10 or \
            choices_parameters['Wgc (µm)']['value'].get() > 10000 or \
            choices_parameters['Lgc (m)']['value'].get() < 0 or \
            choices_parameters['Lgc (m)']['value'].get() > 100:
        messagebox.showerror(title='GC distances', message=
        'GC generally have a thickness and a width between 10µm and 10mm.\
        Also, GC length is generally between 0 and 100m')
        choices.clear()
        return
    if choices_parameters['ε_gdl']['value'].get() < 0 or choices_parameters['ε_gdl']['value'].get() > 1 or \
            choices_parameters['ε_mc']['value'].get() < 0 or choices_parameters['ε_mc']['value'].get() > 1:
        messagebox.showerror(title='Porosities', message='All porosities should be between 0 and 1.')
        choices.clear()
        return
    if choices_parameters['τ']['value'].get() < 1 or choices_parameters['τ']['value'].get() > 4:
        messagebox.showerror(title='Pore structure coefficient', message=
        'The pore structure coefficient should be between 1 and 4.')
        choices.clear()
        return
    if choices_parameters['ε_c']['value'].get() < 0 or choices_parameters['ε_c']['value'].get() > 1:
        messagebox.showerror(title='Compression ratio', message='The compression ratio should be between 0 and 1.')
        choices.clear()
        return
    if choices_parameters['e']['value'].get() < 3 or choices_parameters['e']['value'].get() > 5:
        messagebox.showerror(title='Capillary exponent', message=
        'The capillary exponent should be between 3 and 5 and being an integer.')
        choices.clear()
        return
    if choices_parameters['Re (µΩ.m²)']['value'].get() < 0.5 or choices_parameters['Re (µΩ.m²)']['value'].get() > 5:
        messagebox.showerror(title='Electron conduction resistance', message=
        'The electron conduction resistance is generally between 0.5 and 5 µΩ.m².')
        choices.clear()
        return
    if choices_parameters['i0_c_ref (A/m²)']['value'].get() < 0.001 or \
            choices_parameters['i0_c_ref (A/m²)']['value'].get() > 500:
        messagebox.showerror(title='Referenced exchange current density', message=
        'The referenced exchange current density is generally between 0.001 and 500 A.m-2.')
        choices.clear()
        return
    if choices_parameters['κ_co (mol/(m.s.Pa))']['value'].get() < 0.01 or \
            choices_parameters['κ_co (mol/(m.s.Pa))']['value'].get() > 100:
        messagebox.showerror(title='Crossover correction coefficient', message=
        'The crossover correction coefficient is generally between 0.01 and 100 mol.m-1.s-1.Pa-1.')
        choices.clear()
        return
    if choices_parameters['κ_c']['value'].get() < 0 or choices_parameters['κ_c']['value'].get() > 100:
        messagebox.showerror(title='Overpotential correction exponent', message=
        'The overpotential correction exponent is generally between 0 and 100.')
        choices.clear()
        return
    if choices_parameters['a_slim']['value'].get() < 0 or choices_parameters['a_slim']['value'].get() > 1:
        messagebox.showerror(title='Slop of slim function', message=
        'The slop of slim function is generally between 0 and 1.')
        choices.clear()
        return
    if choices_parameters['b_slim']['value'].get() < 0 or choices_parameters['b_slim']['value'].get() > 1:
        messagebox.showerror(title='Intercept of slim function', message=
        'The intercept of slim function is generally between 0 and 1.')
        choices.clear()
        return
    if choices_parameters['a_switch']['value'].get() < 0 or choices_parameters['a_switch']['value'].get() > 1:
        messagebox.showerror(title='Slop of switch function', message=
        'The slop of switch function is generally between 0 and 1.')
        choices.clear()
        return
    if choices_parameters['C_dl (MF/m³)']['value'].get() < 5 or choices_parameters['C_dl (MF/m³)']['value'].get() > 100:
        messagebox.showerror(title='Double layer capacitance', message='I have not settled yet a range for C_dl.')
        choices.clear()
        return
    if (choices_parameters['t0_step (s)']['value'].get() < 0 or choices_parameters['tf_step (s)']['value'].get() < 0 or \
            choices_parameters['Δt_load_step (s)']['value'].get() < 0 or \
            choices_parameters['Δt_dyn_step (s)']['value'].get() < 0 or \
            choices_parameters['Δt_load_pola (s)']['value'].get() < 0 or \
            choices_parameters['Δt_break_pola (s)']['value'].get() < 0 or \
            choices_parameters['Δt_ini_pola (s)']['value'].get() < 0 or \
            choices_parameters['t0_step (s)']['value'].get() > choices_parameters['tf_step (s)']['value'].get() or \
            choices_parameters['Δt_load_step (s)']['value'].get() >
          (choices_parameters['tf_step (s)']['value'].get() - choices_parameters['t0_step (s)']['value'].get())):
        messagebox.showerror(title='Times', message=
        'The times should be positive, t0_step < tf_step and delta_t_load_step < (tf_step - t0_step).')
        choices.clear()
        return
    if choices_parameters['i_ini_step (A/cm²)']['value'].get() < 0 or \
            choices_parameters['i_final_step (A/cm²)']['value'].get() < 0 or \
            choices_parameters['i_max_pola (A/cm²)']['value'].get() < 0 or \
            choices_parameters['Δi_pola (A/cm²)']['value'].get() < 0 or \
            choices_parameters['i_EIS (A/cm²)']['value'].get() < 0 or \
            choices_parameters['Δi_pola (A/cm²)']['value'].get() > choices_parameters['i_max_pola (A/cm²)']['value'].get() or \
            choices_parameters['i_ini_step (A/cm²)']['value'].get() > choices_parameters['i_final_step (A/cm²)']['value'].get():
        messagebox.showerror(title='Current densities', message=
        'The current densities should be positive, delta_i_pola < i_max_pola and i_ini_step < i_final_step.')
        choices.clear()
        return
    if choices_parameters['ratio_EIS (%)']['value'].get() < 0 or \
            choices_parameters['ratio_EIS (%)']['value'].get() > 20:
        messagebox.showerror(title='Ratio EIS', message=
        'Ratio EIS is a pourcentage of i_EIS and should be between 0 and 20 for plotting correct EIS.')
        choices.clear()
        return

    if choices_parameters['nb_f_EIS']['value'].get() < 0 or \
            choices_parameters['nb_points_EIS']['value'].get() < 0 or \
            type(choices_parameters['f_power_min_EIS']['value'].get()) != int or \
            type(choices_parameters['f_power_max_EIS']['value'].get()) != int or \
            type(choices_parameters['nb_f_EIS']['value'].get()) != int or \
            type(choices_parameters['nb_points_EIS']['value'].get()) != int:
        messagebox.showerror(title='f EIS', message=
        'f_EIS parameters should be integer and number of points should be positive.')
        choices.clear()
        return

    if choices_parameters['t_purge (s)']['value'].get() < 0 or choices_parameters['Δt_purge (s)']['value'].get() < 0:
        messagebox.showerror(title='Purge times', message='Negative times does not characterise purges.')
        choices.clear()
        return

    if choices_parameters['max_step (s)']['value'].get() < 0 or choices_parameters['max_step (s)']['value'].get() > 0.1:
        messagebox.showerror(title='Max step', message=
        'The max step value for the solver should be positive and lower than 0.1 for normal use.')
        choices.clear()
        return

    if choices_parameters['n_gdl']['value'].get() < 5 or type(choices_parameters['n_gdl']['value'].get()) != int:
        messagebox.showerror(title='n gdl', message='The n_gdl value should be an integer bigger than 5.')
        choices.clear()
        return

    if current_button == 0 and choices_buttons['type_display']['value'].get() == 2:
        messagebox.showerror(title='n gdl', message=
        'dynamic plot is not thought to be used with step current and multiple display.' +
        'There would be too much plots to handle.')
        choices.clear()
        return

    if current_button == 2 and choices_buttons['type_plot']['value'].get() == 0:
        messagebox.showerror(title='n gdl', message='EIS has to be plot with a dynamic type_plot setting, '
        'because max_step has to be adjusted at each frequency.')
        choices.clear()
        return

def launch_AlphaPEM_for_step_current(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                                     i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc,
                                     Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co,
                                     kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell,
                                     type_current, type_auxiliary, type_control, type_purge, type_display, type_plot):
    """Launch the AlphaPEM simulator for a step current density and display the results.

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
    """

    # Figures preparation
    fig1, ax1, fig2, ax2 = figures_preparation(type_current, type_display)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if type_plot == "dynamic":
        # Initialization
        #       ... of the plot update number (n) and the initial time interval (time_interval)
        initial_variable_values = None
        t0_step, tf_step, delta_t_load_step, delta_t_dyn_step = t_step
        n = int(tf_step / delta_t_dyn_step)  # It is the plot update number.
        time_interval = [0, delta_t_dyn_step]  # It is the initial time interval.

        # Dynamic simulation
        for i in range(n):
            Simulator = AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                                 i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc,
                                 Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c,
                                 a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell, type_current,
                                 type_auxiliary, type_control, type_purge, type_display, type_plot,
                                 initial_variable_values, time_interval)

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
            if type_display != "no_display":
                Simulator.Display(ax1, ax2)

    else:  # elif type_plot == "fixed":
        # Simulation
        Simulator = AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                             i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                             Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                             b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell, type_current,
                             type_auxiliary, type_control, type_purge, type_display, type_plot)
        # Display
        if type_display != "no_display":
            Simulator.Display(ax1, ax2)
        # Plot saving
        plot_saving(type_fuel_cell, type_current, type_display, fig1, fig2)



def launch_AlphaPEM_for_polarization_current(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step,
                                             i_step, i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                             Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re,
                                             i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step,
                                             n_gdl, t_purge, type_fuel_cell, type_current, type_auxiliary, type_control,
                                             type_purge, type_display, type_plot):
    """Launch the AlphaPEM simulator for a polarization current density and display the results.

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
    """

    # Figures preparation
    fig1, ax1, fig2, ax2 = figures_preparation(type_current, type_display)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if type_plot == "dynamic":
        # Initialization
        #       ... of the plot update number (n) and the initial time interval (time_interval)
        initial_variable_values = None
        delta_t_load_pola, delta_t_break_pola, delta_i_pola, delta_t_ini_pola = delta_pola
        delta_t = delta_t_load_pola + delta_t_break_pola  # s. It is the time of one load.
        tf = delta_t_ini_pola + int(i_max_pola_1 / delta_i_pola + 1) * delta_t  # s. It is the polarization current
        #                                                                            duration.
        n = int(tf / delta_t)  # It is the plot update number.
        time_interval = [0, delta_t_ini_pola + delta_t]  # It is the initial time interval.

        # Dynamic simulation
        for i in range(n):
            Simulator = AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                                 i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc,
                                 Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c,
                                 a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell, type_current,
                                 type_auxiliary, type_control, type_purge, type_display, type_plot,
                                 initial_variable_values, time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                t0_interval = Simulator.variables['t'][-1]
                tf_interval = delta_t_ini_pola + (i + 2) * delta_t
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator.solver_variable_names:
                initial_variable_values.append(Simulator.variables[x][-1])

            # Display
            if type_display != "no_display":
                Simulator.Display(ax1, ax2)

    else:  # elif type_plot == "fixed":
        # Simulation
        Simulator = AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                             i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                             Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                             b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell, type_current,
                             type_auxiliary, type_control, type_purge, type_display, type_plot)
        # Display
        if type_display != "no_display":
            Simulator.Display(ax1, ax2)
        # Plot saving
        plot_saving(type_fuel_cell, type_current, type_display, fig1, fig2)

def launch_AlphaPEM_for_EIS_current(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                                    i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc,
                                    Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co,
                                    kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell,
                                    type_current, type_auxiliary, type_control, type_purge, type_display, type_plot):
    """Launch the AlphaPEM simulator for an EIS current density and display the results.

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
    """

    # Figures preparation
    fig1, ax1, fig2, ax2 = figures_preparation(type_current, type_display)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if type_plot == "dynamic":
        # Initialization
        #       ... of the plot update number (n) and the initial time interval (time_interval)
        initial_variable_values = None
        t0_EIS, t_new_start, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
        f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS  # These are used for EIS max_step
        #                                                                    actualization.
        f = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)  # It is a list of all the frequency tested.
        n = len(t_new_start)  # It is the plot update number.
        time_interval = [0, t0_EIS]  # It is the initial time interval.

        #       A preliminary simulation run is necessary to equilibrate the internal variables of the cell at i_EIS
        #       prior to initiating the EIS.
        Simulator = AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                             i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                             Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                             b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell, type_current,
                             type_auxiliary, type_control, type_purge, type_display, type_plot,
                             initial_variable_values, time_interval)

        # time_interval actualization
        t0_EIS_temp = t0_EIS  # It is the initial time for 1 EIS point.
        tf_EIS_temp = t_new_start[0] + delta_t_break_EIS[0] + delta_t_measurement_EIS[0]  # It is the final time for
        #                                                                                  1 EIS point.
        n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which has been
        #                                                      made.
        max_step = 1 / (f[n_inf] * nb_points_EIS)  # max_step is actualized according to the current frequency
        #                                        for increased calculation
        time_interval = [t0_EIS_temp, tf_EIS_temp]

        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values = []
        for x in Simulator.solver_variable_names:
            initial_variable_values.append(Simulator.variables[x][-1])

        # Dynamic simulation
        for i in range(n):
            Simulator = AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                                 i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc,
                                 Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co,
                                 kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell,
                                 type_current, type_auxiliary, type_control, type_purge, type_display, type_plot,
                                 initial_variable_values, time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                t0_EIS_temp = Simulator.variables['t'][-1]  # It is the initial time for 1 EIS point.
                tf_EIS_temp = t_new_start[i + 1] + delta_t_break_EIS[i + 1] + delta_t_measurement_EIS[i + 1]  # It
                #                                                                 is the final time for 1 EIS point.
                n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which
                #                                                      has been made.
                max_step = 1 / (f[n_inf] * nb_points_EIS)  # max_step is actualized according to the current
                #                                            frequency for increased calculation
                time_interval = [t0_EIS_temp, tf_EIS_temp]  # It is the time interval for 1 EIS point.

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator.solver_variable_names:
                initial_variable_values.append(Simulator.variables[x][-1])

            # Display
            if type_display != "no_display":
                Simulator.Display(ax1, ax2)

    else:  # elif type_plot == "fixed":
        # Simulation
        Simulator = AlphaPEM(current_density, Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, t_step, i_step,
                             i_max_pola, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc,
                             Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co,
                             kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell,
                             type_current, type_auxiliary, type_control, type_purge, type_display, type_plot)
        # Display
        if type_display != "no_display":
            Simulator.Display(ax1, ax2)
        # Plot saving
        plot_saving(type_fuel_cell, type_current, type_display, fig1, fig2)