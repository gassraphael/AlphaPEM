# -*- coding: utf-8 -*-
import numpy as np
import tkinter as tk
from tkinter import messagebox

# Constants value and functions
from modules.settings_modules import stored_operating_inputs, stored_physical_parameters

"""
This modul contains some of the required functions for the main program AlphaPEM_interface.
"""


def changeValue(frame, choices, label_widgets, entry_widgets):
    """
    This function is called when the user selects a specific option from a dropdown menu 
    (choices['setting_input']). Depending on the selected option, it either hides or displays specific 
    input fields (labels or entry widgets) on the GUI.
    """

    if choices['setting_input'].get() == 'EH-31 1.5 bar (2021)' or \
            choices['setting_input'].get() == 'EH-31 2.0 bar (2021)' or \
            choices['setting_input'].get() == 'EH-31 2.25 bar (2021)' or \
            choices['setting_input'].get() == 'EH-31 2.5 bar (2021)' or \
            choices['setting_input'].get() == 'Biao Xie 1.0 bar (2015)' or \
            choices['setting_input'].get() == 'Biao Xie 1.35 bar (2015)' or \
            choices['setting_input'].get() == 'Linhao Fan (2010)':
        # Hide all the Entries if they are already displayed
        for entry in entry_widgets:
            entry.grid_forget()

            # Recovers the new settings
        recover_for_display_operating_inputs_and_physical_parameters(frame, choices)

        # Display them
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Tfc_i']). \
            grid(row=2, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Pa_des_i']). \
            grid(row=2, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Pc_des_i']). \
            grid(row=2, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Sa_i']). \
            grid(row=3, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Sc_i']). \
            grid(row=3, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Phi_a_des_i']). \
            grid(row=4, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Phi_c_des_i']). \
            grid(row=4, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Aact_i']). \
            grid(row=6, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Hgdl_i']). \
            grid(row=7, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Hcl_i']). \
            grid(row=7, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Hmem_i']). \
            grid(row=7, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Hgc_i']). \
            grid(row=8, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Wgc_i']). \
            grid(row=8, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Lgc_i']). \
            grid(row=8, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['epsilon_gdl_i']). \
            grid(row=10, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['epsilon_mc_i']). \
            grid(row=10, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['tau_i']). \
            grid(row=10, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['epsilon_c_i']). \
            grid(row=11, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['e_i']). \
            grid(row=11, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['Re_i']). \
            grid(row=11, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['i0_c_ref_i']). \
            grid(row=12, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['kappa_co_i']). \
            grid(row=12, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['kappa_c_i']). \
            grid(row=12, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['a_slim_i']). \
            grid(row=13, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['b_slim_i']). \
            grid(row=13, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['a_switch_i']). \
            grid(row=13, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['C_dl_i']). \
            grid(row=14, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['t0_step_i']). \
            grid(row=16, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['tf_step_i']). \
            grid(row=16, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['delta_t_load_step_i']). \
            grid(row=16, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['i_ini_step_i']). \
            grid(row=17, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['i_final_step_i']). \
            grid(row=17, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['i_max_pola_i']). \
            grid(row=18, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['delta_i_pola_i']). \
            grid(row=18, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['delta_t_load_pola_i']). \
            grid(row=19, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['delta_t_break_pola_i']). \
            grid(row=19, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['delta_t_ini_pola_i']). \
            grid(row=19, column=5, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['i_EIS_i']). \
            grid(row=20, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['delta_t_dyn_step_i']). \
            grid(row=23, column=1, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w', textvariable=choices['t_purge_i']). \
            grid(row=23, column=3, padx=5)
        tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                 textvariable=choices['delta_t_purge_i']). \
            grid(row=23, column=5, padx=5)

    else:  # setting_input == 'Enter its specifications':
        # Hide the Labels if they are already displayed
        for label in label_widgets:
            label.grid_forget()

        # Saves and displays the user entries
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Tfc_i']). \
            grid(row=2, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Pa_des_i']). \
            grid(row=2, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Pc_des_i']). \
            grid(row=2, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Sa_i']). \
            grid(row=3, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Sc_i']). \
            grid(row=3, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Phi_a_des_i']). \
            grid(row=4, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Phi_c_des_i']). \
            grid(row=4, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Aact_i']). \
            grid(row=6, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Hgdl_i']). \
            grid(row=7, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Hcl_i']). \
            grid(row=7, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Hmem_i']). \
            grid(row=7, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Hgc_i']). \
            grid(row=8, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Wgc_i']). \
            grid(row=8, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Lgc_i']). \
            grid(row=8, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['epsilon_gdl_i']). \
            grid(row=10, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['epsilon_mc_i']). \
            grid(row=10, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['tau_i']). \
            grid(row=10, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['epsilon_c_i']). \
            grid(row=11, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['e_i']). \
            grid(row=11, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['Re_i']). \
            grid(row=11, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['i0_c_ref_i']). \
            grid(row=12, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['kappa_co_i']). \
            grid(row=12, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['kappa_c_i']). \
            grid(row=12, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['a_slim_i']). \
            grid(row=13, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['b_slim_i']). \
            grid(row=13, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['a_switch_i']). \
            grid(row=13, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['C_dl_i']). \
            grid(row=14, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['t0_step_i']). \
            grid(row=16, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['tf_step_i']). \
            grid(row=16, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_load_step_i']). \
            grid(row=16, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['i_ini_step_i']). \
            grid(row=17, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['i_final_step_i']). \
            grid(row=17, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['i_max_pola_i']). \
            grid(row=18, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_i_pola_i']). \
            grid(row=18, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_load_pola_i']). \
            grid(row=19, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_break_pola_i']). \
            grid(row=19, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_ini_pola_i']). \
            grid(row=19, column=5, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['i_EIS_i']). \
            grid(row=20, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_dyn_step_i']). \
            grid(row=23, column=1, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['t_purge_i']). \
            grid(row=23, column=3, padx=5)
        tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_purge_i']). \
            grid(row=23, column=5, padx=5)


def display_label_operating_inputs_and_physical_parameters(frame):
    """
    This function displays labels on the GUI, representing operating conditions and physical parameters, 
    without their actual values.
    """

    # Operating conditions
    tk.Label(frame, text='Operating conditions', fg='black', font=('Times New Roman', 12, 'bold')). \
        grid(row=1, column=0, columnspan=6, ipady=15)
    tk.Label(frame, text='Tfc (°C)', fg='black', font=('Times New Roman', 10)).grid(row=2, column=0, sticky="w")
    tk.Label(frame, text='Pa_des (bar)', fg='black', font=('Times New Roman', 10)).grid(row=2, column=2, sticky="w")
    tk.Label(frame, text='Pc_des (bar)', fg='black', font=('Times New Roman', 10)).grid(row=2, column=4, sticky="w")
    tk.Label(frame, text='Sa', fg='black', font=('Times New Roman', 10)).grid(row=3, column=0, sticky="w")
    tk.Label(frame, text='Sc', fg='black', font=('Times New Roman', 10)).grid(row=3, column=2, sticky="w")
    tk.Label(frame, text='Ф_a_des', fg='black', font=('Times New Roman', 10)).grid(row=4, column=0, sticky="w")
    tk.Label(frame, text='Ф_c_des', fg='black', font=('Times New Roman', 10)).grid(row=4, column=2, sticky="w")

    # Accessible physical parameters
    tk.Label(frame, text='Accessible physical parameters', fg='black', font=('Times New Roman', 12, 'bold')). \
        grid(row=5, column=0, columnspan=6, ipady=15)
    tk.Label(frame, text='Aact (cm²)', fg='black', font=('Times New Roman', 10)).grid(row=6, column=0, sticky="w")
    tk.Label(frame, text='Hgdl (µm)', fg='black', font=('Times New Roman', 10)).grid(row=7, column=0, sticky="w")
    tk.Label(frame, text='Hcl (µm)', fg='black', font=('Times New Roman', 10)).grid(row=7, column=2, sticky="w")
    tk.Label(frame, text='Hmem (µm)', fg='black', font=('Times New Roman', 10)).grid(row=7, column=4, sticky="w")
    tk.Label(frame, text='Hgc (µm)', fg='black', font=('Times New Roman', 10)).grid(row=8, column=0, sticky="w")
    tk.Label(frame, text='Wgc (µm)', fg='black', font=('Times New Roman', 10)).grid(row=8, column=2, sticky="w")
    tk.Label(frame, text='Lgc (m)', fg='black', font=('Times New Roman', 10)).grid(row=8, column=4, sticky="w")

    # Undetermined physical parameters
    tk.Label(frame, text='Undetermined physical parameters', fg='black', font=('Times New Roman', 12, 'bold')). \
        grid(row=9, column=0, columnspan=6, ipady=15)
    tk.Label(frame, text='ε_gdl', fg='black', font=('Times New Roman', 10)).grid(row=10, column=0, sticky="w")
    tk.Label(frame, text='ε_mc', fg='black', font=('Times New Roman', 10)).grid(row=10, column=2, sticky="w")
    tk.Label(frame, text='τ', fg='black', font=('Times New Roman', 10)).grid(row=10, column=4, sticky="w")
    tk.Label(frame, text='ε_c', fg='black', font=('Times New Roman', 10)).grid(row=11, column=0, sticky="w")
    tk.Label(frame, text='e', fg='black', font=('Times New Roman', 10)).grid(row=11, column=2, sticky="w")
    tk.Label(frame, text='Re (µΩ.m²)', fg='black', font=('Times New Roman', 10)).grid(row=11, column=4, sticky="w")
    tk.Label(frame, text='i0_c_ref (A/m²)', fg='black', font=('Times New Roman', 10)). \
        grid(row=12, column=0, sticky="w")
    tk.Label(frame, text='κ_co (mol/(m.s.Pa))', fg='black', font=('Times New Roman', 10)). \
        grid(row=12, column=2, sticky="w")
    tk.Label(frame, text='κ_c', fg='black', font=('Times New Roman', 10)).grid(row=12, column=4, sticky="w")
    tk.Label(frame, text='a_slim', fg='black', font=('Times New Roman', 10)).grid(row=13, column=0, sticky="w")
    tk.Label(frame, text='b_slim', fg='black', font=('Times New Roman', 10)).grid(row=13, column=2, sticky="w")
    tk.Label(frame, text='a_switch', fg='black', font=('Times New Roman', 10)).grid(row=13, column=4, sticky="w")
    tk.Label(frame, text='C_dl (MF/m³)', fg='black', font=('Times New Roman', 10)).grid(row=14, column=0, sticky="w")

    # Current density parameters
    tk.Label(frame, text='Current density parameters', fg='black', font=('Times New Roman', 12, 'bold')). \
        grid(row=15, column=0, columnspan=6, ipady=15)
    tk.Label(frame, text='t0_step (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=16, column=0, sticky="w")
    tk.Label(frame, text='tf_step (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=16, column=2, sticky="w")
    tk.Label(frame, text='Δt_load_step (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=16, column=4, sticky="w")
    tk.Label(frame, text='i_ini_step (A/cm²)', fg='black', font=('Times New Roman', 10)). \
        grid(row=17, column=0, sticky="w")
    tk.Label(frame, text='i_final_step (A/cm²)', fg='black', font=('Times New Roman', 10)). \
        grid(row=17, column=2, sticky="w")
    tk.Label(frame, text='i_max_pola (A/cm²)', fg='black', font=('Times New Roman', 10)). \
        grid(row=18, column=0, sticky="w")
    tk.Label(frame, text='Δi_pola (A/cm²)', fg='black', font=('Times New Roman', 10)). \
        grid(row=18, column=2, sticky="w")
    tk.Label(frame, text='Δt_load_pola (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=19, column=0, sticky="w")
    tk.Label(frame, text='Δt_break_pola (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=19, column=2, sticky="w")
    tk.Label(frame, text='Δt_ini_pola (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=19, column=4, sticky="w")
    tk.Label(frame, text='i_EIS (A/cm²)', fg='black', font=('Times New Roman', 10)). \
        grid(row=20, column=0, sticky="w")

    # Computing parameters
    tk.Label(frame, text='Computing parameters', fg='black', font=('Times New Roman', 12, 'bold')). \
        grid(row=22, column=0, columnspan=6, ipady=15)
    tk.Label(frame, text='Δt_dyn_step (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=23, column=0, sticky="w")
    tk.Label(frame, text='t_purge (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=23, column=2, sticky="w")
    tk.Label(frame, text='Δt_purge (s)', fg='black', font=('Times New Roman', 10)). \
        grid(row=23, column=4, sticky="w")


def display_value_operating_inputs_and_physical_parameters(frame, choices):
    """
    This function displays entry widgets on the GUI, 
    where the user can enter values for operating conditions and physical parameters.
    """

    Label_Tfc = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                         textvariable=choices['Tfc_i'])
    Label_Pa_des = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                            textvariable=choices['Pa_des_i'])
    Label_Pc_des = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                            textvariable=choices['Pc_des_i'])
    Label_Sa = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                        textvariable=choices['Sa_i'])
    Label_Sc = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                        textvariable=choices['Sc_i'])
    Label_Phi_a_des = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                               textvariable=choices['Phi_a_des_i'])
    Label_Phi_c_des = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                               textvariable=choices['Phi_c_des_i'])
    Label_Aact = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                          textvariable=choices['Aact_i'])
    Label_Hgdl = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                          textvariable=choices['Hgdl_i'])
    Label_Hcl = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                         textvariable=choices['Hcl_i'])
    Label_Hmem = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                          textvariable=choices['Hmem_i'])
    Label_Hgc = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                         textvariable=choices['Hgc_i'])
    Label_Wgc = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                         textvariable=choices['Wgc_i'])
    Label_Lgc = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                         textvariable=choices['Lgc_i'])
    Label_epsilon_gdl = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                 textvariable=choices['epsilon_gdl_i'])
    Label_epsilon_mc = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                textvariable=choices['epsilon_mc_i'])
    Label_tau = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                         textvariable=choices['tau_i'])
    Label_epsilon_c = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                               textvariable=choices['epsilon_c_i'])
    Label_e = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                       textvariable=choices['e_i'])
    Label_Re = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                        textvariable=choices['Re_i'])
    Label_i0_c_ref = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                              textvariable=choices['i0_c_ref_i'])
    Label_kappa_co = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                              textvariable=choices['kappa_co_i'])
    Label_kappa_c = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                             textvariable=choices['kappa_c_i'])
    Label_a_slim = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                            textvariable=choices['a_slim_i'])
    Label_b_slim = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                            textvariable=choices['b_slim_i'])
    Label_a_switch = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                              textvariable=choices['a_switch_i'])
    Label_C_dl = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                          textvariable=choices['C_dl_i'])
    Label_t0_step = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                             textvariable=choices['t0_step_i'])
    Label_tf_step = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                             textvariable=choices['tf_step_i'])
    Label_delta_t_load_step = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                       textvariable=choices['delta_t_load_step_i'])
    Label_i_ini_step = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                textvariable=choices['i_ini_step_i'])
    Label_i_final_step = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                  textvariable=choices['i_final_step_i'])
    Label_i_max_pola = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                textvariable=choices['i_max_pola_i'])
    Label_delta_i_pola = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                  textvariable=choices['delta_i_pola_i'])
    Label_delta_t_load_pola = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                       textvariable=choices['delta_t_load_pola_i'])
    Label_delta_t_break_pola = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                        textvariable=choices['delta_t_break_pola_i'])
    Label_delta_t_ini_pola = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                      textvariable=choices['delta_t_ini_pola_i'])
    Label_i_EIS = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                           textvariable=choices['i_EIS_i'])
    Label_delta_t_dyn_step = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                      textvariable=choices['delta_t_dyn_step_i'])
    Label_t_purge = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                             textvariable=choices['t_purge_i'])
    Label_delta_t_purge = tk.Label(frame, width=7, bg='white', bd=1, relief='solid', anchor='w',
                                   textvariable=choices['delta_t_purge_i'])

    label_widgets = [Label_Tfc, Label_Pa_des, Label_Pc_des, Label_Sa, Label_Sc, Label_Phi_a_des, Label_Phi_c_des,
                     Label_Aact, Label_Hgdl, Label_Hcl, Label_Hmem, Label_Hgc, Label_Wgc, Label_Lgc, Label_epsilon_gdl,
                     Label_epsilon_mc, Label_tau, Label_epsilon_c, Label_e, Label_Re, Label_i0_c_ref, Label_kappa_co,
                     Label_kappa_c, Label_a_slim, Label_b_slim, Label_a_switch, Label_C_dl, Label_t0_step,
                     Label_tf_step, Label_delta_t_load_step, Label_i_ini_step, Label_i_final_step, Label_i_max_pola,
                     Label_delta_i_pola, Label_delta_t_load_pola, Label_delta_t_break_pola, Label_delta_t_ini_pola,
                     Label_i_EIS, Label_delta_t_dyn_step, Label_t_purge, Label_delta_t_purge]

    Entry_Tfc = tk.Entry(frame, width=7, bg='white', textvariable=choices['Tfc_i'])
    Entry_Pa_des = tk.Entry(frame, width=7, bg='white', textvariable=choices['Pa_des_i'])
    Entry_Pc_des = tk.Entry(frame, width=7, bg='white', textvariable=choices['Pc_des_i'])
    Entry_Sa = tk.Entry(frame, width=7, bg='white', textvariable=choices['Sa_i'])
    Entry_Sc = tk.Entry(frame, width=7, bg='white', textvariable=choices['Sc_i'])
    Entry_Phi_a_des = tk.Entry(frame, width=7, bg='white', textvariable=choices['Phi_a_des_i'])
    Entry_Phi_c_des = tk.Entry(frame, width=7, bg='white', textvariable=choices['Phi_c_des_i'])
    Entry_Aact = tk.Entry(frame, width=7, bg='white', textvariable=choices['Aact_i'])
    Entry_Hgdl = tk.Entry(frame, width=7, bg='white', textvariable=choices['Hgdl_i'])
    Entry_Hcl = tk.Entry(frame, width=7, bg='white', textvariable=choices['Hcl_i'])
    Entry_Hmem = tk.Entry(frame, width=7, bg='white', textvariable=choices['Hmem_i'])
    Entry_Hgc = tk.Entry(frame, width=7, bg='white', textvariable=choices['Hgc_i'])
    Entry_Wgc = tk.Entry(frame, width=7, bg='white', textvariable=choices['Wgc_i'])
    Entry_Lgc = tk.Entry(frame, width=7, bg='white', textvariable=choices['Lgc_i'])
    Entry_epsilon_gdl = tk.Entry(frame, width=7, bg='white', textvariable=choices['epsilon_gdl_i'])
    Entry_epsilon_mc = tk.Entry(frame, width=7, bg='white', textvariable=choices['epsilon_mc_i'])
    Entry_tau = tk.Entry(frame, width=7, bg='white', textvariable=choices['tau_i'])
    Entry_epsilon_c = tk.Entry(frame, width=7, bg='white', textvariable=choices['epsilon_c_i'])
    Entry_e = tk.Entry(frame, width=7, bg='white', textvariable=choices['e_i'])
    Entry_Re = tk.Entry(frame, width=7, bg='white', textvariable=choices['Re_i'])
    Entry_i0_c_ref = tk.Entry(frame, width=7, bg='white', textvariable=choices['i0_c_ref_i'])
    Entry_kappa_co = tk.Entry(frame, width=7, bg='white', textvariable=choices['kappa_co_i'])
    Entry_kappa_c = tk.Entry(frame, width=7, bg='white', textvariable=choices['kappa_c_i'])
    Entry_a_slim = tk.Entry(frame, width=7, bg='white', textvariable=choices['a_slim_i'])
    Entry_b_slim = tk.Entry(frame, width=7, bg='white', textvariable=choices['b_slim_i'])
    Entry_a_switch = tk.Entry(frame, width=7, bg='white', textvariable=choices['a_switch_i'])
    Entry_C_dl = tk.Entry(frame, width=7, bg='white', textvariable=choices['C_dl_i'])
    Entry_t0_step = tk.Entry(frame, width=7, bg='white', textvariable=choices['t0_step_i'])
    Entry_tf_step = tk.Entry(frame, width=7, bg='white', textvariable=choices['tf_step_i'])
    Entry_delta_t_load_step = tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_load_step_i'])
    Entry_i_ini_step = tk.Entry(frame, width=7, bg='white', textvariable=choices['i_ini_step_i'])
    Entry_i_final_step = tk.Entry(frame, width=7, bg='white', textvariable=choices['i_final_step_i'])
    Entry_i_max_pola = tk.Entry(frame, width=7, bg='white', textvariable=choices['i_max_pola_i'])
    Entry_delta_i_pola = tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_i_pola_i'])
    Entry_delta_t_load_pola = tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_load_pola_i'])
    Entry_delta_t_break_pola = tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_break_pola_i'])
    Entry_delta_t_ini_pola = tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_ini_pola_i'])
    Entry_i_EIS = tk.Entry(frame, width=7, bg='white', textvariable=choices['i_EIS_i'])
    Entry_delta_t_dyn_step = tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_dyn_step_i'])
    Entry_t_purge = tk.Entry(frame, width=7, bg='white', textvariable=choices['t_purge_i'])
    Entry_delta_t_purge = tk.Entry(frame, width=7, bg='white', textvariable=choices['delta_t_purge_i'])

    entry_widgets = [Entry_Tfc, Entry_Pa_des, Entry_Pc_des, Entry_Sa, Entry_Sc, Entry_Phi_a_des, Entry_Phi_c_des,
                     Entry_Aact, Entry_Hgdl, Entry_Hcl, Entry_Hmem, Entry_Hgc, Entry_Wgc, Entry_Lgc, Entry_epsilon_gdl,
                     Entry_epsilon_mc, Entry_tau, Entry_epsilon_c, Entry_e, Entry_Re, Entry_i0_c_ref, Entry_kappa_co,
                     Entry_kappa_c, Entry_a_slim, Entry_b_slim, Entry_a_switch, Entry_C_dl, Entry_t0_step,
                     Entry_tf_step, Entry_delta_t_load_step, Entry_i_ini_step, Entry_i_final_step, Entry_i_max_pola,
                     Entry_i_max_pola, Entry_delta_t_load_pola, Entry_delta_t_break_pola, Entry_delta_t_ini_pola,
                     Entry_i_EIS, Entry_delta_t_dyn_step, Entry_t_purge, Entry_delta_t_purge]

    Entry_Tfc.grid(row=2, column=1, padx=5)
    Entry_Pa_des.grid(row=2, column=3, padx=5)
    Entry_Pc_des.grid(row=2, column=5, padx=5)
    Entry_Sa.grid(row=3, column=1, padx=5)
    Entry_Sc.grid(row=3, column=3, padx=5)
    Entry_Phi_a_des.grid(row=4, column=1, padx=5)
    Entry_Phi_c_des.grid(row=4, column=3, padx=5)
    Entry_Aact.grid(row=6, column=1, padx=5)
    Entry_Hgdl.grid(row=7, column=1, padx=5)
    Entry_Hcl.grid(row=7, column=3, padx=5)
    Entry_Hmem.grid(row=7, column=5, padx=5)
    Entry_Hgc.grid(row=8, column=1, padx=5)
    Entry_Wgc.grid(row=8, column=3, padx=5)
    Entry_Lgc.grid(row=8, column=5, padx=5)
    Entry_epsilon_gdl.grid(row=10, column=1, padx=5)
    Entry_epsilon_mc.grid(row=10, column=3, padx=5)
    Entry_tau.grid(row=10, column=5, padx=5)
    Entry_epsilon_c.grid(row=11, column=1, padx=5)
    Entry_e.grid(row=11, column=3, padx=5)
    Entry_Re.grid(row=11, column=5, padx=5)
    Entry_i0_c_ref.grid(row=12, column=1, padx=5)
    Entry_kappa_co.grid(row=12, column=3, padx=5)
    Entry_kappa_c.grid(row=12, column=5, padx=5)
    Entry_a_slim.grid(row=13, column=1, padx=5)
    Entry_b_slim.grid(row=13, column=3, padx=5)
    Entry_a_switch.grid(row=13, column=5, padx=5)
    Entry_C_dl.grid(row=14, column=1, padx=5)
    Entry_t0_step.grid(row=16, column=1, padx=5)
    Entry_tf_step.grid(row=16, column=3, padx=5)
    Entry_delta_t_load_step.grid(row=16, column=5, padx=5)
    Entry_i_ini_step.grid(row=17, column=1, padx=5)
    Entry_i_final_step.grid(row=17, column=3, padx=5)
    Entry_i_max_pola.grid(row=18, column=1, padx=5)
    Entry_delta_i_pola.grid(row=18, column=3, padx=5)
    Entry_delta_t_load_pola.grid(row=19, column=1, padx=5)
    Entry_delta_t_break_pola.grid(row=19, column=3, padx=5)
    Entry_delta_t_ini_pola.grid(row=19, column=5, padx=5)
    Entry_i_EIS.grid(row=20, column=1, padx=5)
    Entry_delta_t_dyn_step.grid(row=23, column=1, padx=5)
    Entry_t_purge.grid(row=23, column=3, padx=5)
    Entry_delta_t_purge.grid(row=23, column=5, padx=5)

    return label_widgets, entry_widgets


def display_radiobuttons(frame, choices):
    """
    This function displays radiobuttons on the GUI, allowing the user to make choices for purging, 
    results display, plot style, etc.
    """
    tk.Label(frame, text='Model possibilities', fg='black', font=('Times New Roman', 12, 'bold')) \
        .grid(row=24, column=0, columnspan=6, ipady=15)

    # Ask the user to choose an option and save it
    tk.Label(frame, text='Auxiliaries: ', fg='black', font=('Times New Roman', 12)). \
        grid(row=25, column=0, sticky="w")
    tk.Radiobutton(frame, text='No auxiliaries', value=2, variable=choices['auxiliaries_choice']). \
        grid(row=25, column=1, sticky="w")
    tk.Radiobutton(frame, text='Closed anode', value=1, variable=choices['auxiliaries_choice']). \
        grid(row=25, column=2, sticky="w")
    tk.Radiobutton(frame, text='Opened anode', value=0, variable=choices['auxiliaries_choice']). \
        grid(row=25, column=3, sticky="w")

    # Ask the user to choose an option and save it
    tk.Label(frame, text='Purge: ', fg='black', font=('Times New Roman', 12)).grid(row=26, column=0, sticky="w")
    tk.Radiobutton(frame, text='Periodic', value=1, variable=choices['is_purging']).grid(row=26, column=1, sticky="w")
    tk.Radiobutton(frame, text='No purge', value=0, variable=choices['is_purging']).grid(row=26, column=2, sticky="w")

    # Ask the user to choose an option and save it
    tk.Label(frame, text='Results: ', fg='black', font=('Times New Roman', 12)).grid(row=27, column=0, sticky="w")
    tk.Radiobutton(frame, text='Precise', value=1, variable=choices['is_precise']).grid(row=27, column=1, sticky="w")
    tk.Radiobutton(frame, text='Fast', value=0, variable=choices['is_precise']).grid(row=27, column=2, sticky="w")

    # Ask the user to choose an option and save it
    tk.Label(frame, text='Display: ', fg='black', font=('Times New Roman', 12)). \
        grid(row=28, column=0, sticky="w")
    tk.Radiobutton(frame, text='Synthetic', value=1, variable=choices['is_synthetic']). \
        grid(row=28, column=1, sticky="w")
    tk.Radiobutton(frame, text='Shattered', value=0, variable=choices['is_synthetic']). \
        grid(row=28, column=2, sticky="w")

    # Ask the user to choose an option and save it
    tk.Label(frame, text='Plot: ', fg='black', font=('Times New Roman', 12)).grid(row=29, column=0, sticky="w")
    tk.Radiobutton(frame, text='Dynamic', value=1, variable=choices['is_dynamic']).grid(row=29, column=1, sticky="w")
    tk.Radiobutton(frame, text='Fixed', value=0, variable=choices['is_dynamic']).grid(row=29, column=2, sticky="w")


def recover_for_display_operating_inputs_and_physical_parameters(frame, choices):
    """
    This function retrieves parameter values for predefined stacks (e.g., "EH-31 1.5 bar", "Biao Xie 1.0 bar", etc.)
    and converts them to appropriate units for display on the GUI.
    """

    if choices['setting_input'].get() == "EH-31 1.5 bar (2021)":
        type_fuel_cell = "EH-31_1.5"
    elif choices['setting_input'].get() == "EH-31 2.0 bar (2021)":
        type_fuel_cell = "EH-31_2.0"
    elif choices['setting_input'].get() == "EH-31 2.25 bar (2021)":
        type_fuel_cell = "EH-31_2.25"
    elif choices['setting_input'].get() == "EH-31 2.5 bar (2021)":
        type_fuel_cell = "EH-31_2.5"
    elif choices['setting_input'].get() == "Biao Xie 1.0 bar (2015)":
        type_fuel_cell = "BX_1.0"
    elif choices['setting_input'].get() == "Biao Xie 1.35 bar (2015)":
        type_fuel_cell = "BX_1.35"
    elif choices['setting_input'].get() == "Linhao Fan (2010)":
        type_fuel_cell = "LF"
    else:
        raise ValueError('the type_fuel_cell given is not valid.')

    Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, i_max_pola = stored_operating_inputs(type_fuel_cell)

    Hcl, epsilon_mc, tau, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, Re, i0_c_ref, kappa_co, kappa_c, \
        a_slim, b_slim, a_switch, C_dl = \
        stored_physical_parameters(type_fuel_cell)

    choices['Tfc_i'].set(np.round(Tfc - 273.15))  # °C
    choices['Pa_des_i'].set(np.round(Pa_des / 1e5, 2))  # bar
    choices['Pc_des_i'].set(np.round(Pc_des / 1e5, 2))  # bar
    choices['Sa_i'].set(np.round(Sa, 1))
    choices['Sc_i'].set(np.round(Sc, 1))
    choices['Phi_a_des_i'].set(np.round(Phi_a_des, 1))
    choices['Phi_c_des_i'].set(np.round(Phi_c_des, 1))
    choices['Aact_i'].set(np.round(Aact * 1e4))  # cm²
    choices['Hgdl_i'].set(np.round(Hgdl * 1e6))  # µm
    choices['Hcl_i'].set(np.round(Hcl * 1e6))  # µm
    choices['Hmem_i'].set(np.round(Hmem * 1e6))  # µm
    choices['Hgc_i'].set(np.round(Hgc * 1e6))  # µm
    choices['Wgc_i'].set(np.round(Wgc * 1e6))  # µm
    choices['Lgc_i'].set(np.round(Lgc, 2))  # m
    choices['epsilon_gdl_i'].set(np.round(epsilon_gdl, 2))
    choices['epsilon_mc_i'].set(np.round(epsilon_mc, 2))
    choices['tau_i'].set(np.round(tau, 2))
    choices['epsilon_c_i'].set(np.round(epsilon_c, 2))
    choices['e_i'].set(e)
    choices['Re_i'].set(np.round(Re * 1e6, 2))  # µΩ.m²
    choices['i0_c_ref_i'].set(np.round(i0_c_ref, 2))  # A.m-2
    choices['kappa_co_i'].set(np.round(kappa_co, 2))  # mol.m-1.s-1.Pa-1
    choices['kappa_c_i'].set(np.round(kappa_c, 2))
    choices['a_slim_i'].set(np.round(a_slim, 4))
    choices['b_slim_i'].set(np.round(b_slim, 4))
    choices['a_switch_i'].set(np.round(a_switch, 4))
    choices['C_dl_i'].set(np.round(C_dl * 1e-6, 2))  # MF.m-3
    choices['i_max_pola_i'].set(np.round(i_max_pola / 1e4, 2))  # A/cm²


def recover_for_use_operating_inputs_and_physical_parameters(choices):
    """
    This function retrieves and converts the parameter values from the GUI into standard units 
    for further calculations.
    """
    Tfc = choices['Tfc_i'].get() + 273.15  # K
    Pa_des, Pc_des = choices['Pa_des_i'].get() * 1e5, choices['Pc_des_i'].get() * 1e5  # Pa
    Sa, Sc = choices['Sa_i'].get(), choices['Sc_i'].get()
    Phi_a_des, Phi_c_des = choices['Phi_a_des_i'].get(), choices['Phi_c_des_i'].get()
    Aact = choices['Aact_i'].get() * 1e-4  # m²
    Hgdl = choices['Hgdl_i'].get() * 1e-6  # m
    Hcl = choices['Hcl_i'].get() * 1e-6  # m
    Hmem = choices['Hmem_i'].get() * 1e-6  # m
    Hgc = choices['Hgc_i'].get() * 1e-6  # m
    Wgc = choices['Wgc_i'].get() * 1e-6  # m
    Lgc = choices['Lgc_i'].get()  # m
    epsilon_gdl = choices['epsilon_gdl_i'].get()
    epsilon_mc = choices['epsilon_mc_i'].get()
    tau = choices['tau_i'].get()
    epsilon_c = choices['epsilon_c_i'].get()
    e = choices['e_i'].get()
    Re = choices['Re_i'].get() * 1e-6  # ohm.m²
    i0_c_ref = choices['i0_c_ref_i'].get()  # A.m-2
    kappa_co = choices['kappa_co_i'].get()  # mol.m-1.s-1.Pa-1
    kappa_c = choices['kappa_c_i'].get()
    a_slim = choices['a_slim_i'].get()
    b_slim = choices['b_slim_i'].get()
    a_switch = choices['a_switch_i'].get()
    C_dl = choices['C_dl_i'].get() * 1e6  # F.m-3
    t_step = (choices['t0_step_i'].get(), choices['tf_step_i'].get(),
              choices['delta_t_load_step_i'].get(), choices['delta_t_dyn_step_i'].get())  # (s, s, s, s)
    i_step = choices['i_ini_step_i'].get() * 1e4, choices['i_final_step_i'].get() * 1e4  # (A.m-2, A.m-2)
    i_max_pola = choices['i_max_pola_i'].get() * 1e4  # A.m-2
    delta_pola = choices['delta_t_load_pola'].get(), choices['delta_t_break_pola'].get(), \
        choices['delta_i_pola'].get() * 1e4, choices['delta_t_ini_pola'].get()  # (s, s, A.m-2, s)

    i_EIS = choices['i_EIS_i'].get() * 1e4

    t_purge, delta_t_purge = choices['t_purge_i'].get(), choices['delta_t_purge_i'].get()  # s

    if choices['setting_input'].get() == "EH-31 1.5 bar":
        type_fuel_cell = "EH-31_1.5"
    elif choices['setting_input'].get() == "Biao Xie 1.0 bar":
        type_fuel_cell = "BX_1.0"
    elif choices['setting_input'].get() == "Linhao Fan":
        type_fuel_cell = "LF"
    elif choices['setting_input'].get() == "Manual setup":
        type_fuel_cell = "manual_setup"

    if choices['auxiliaries_choice'].get() == 0:
        type_auxiliary = "opened_anode"
    elif choices['auxiliaries_choice'].get() == 1:
        type_auxiliary = "closed_anode"
    else:
        type_auxiliary = "no_auxiliary"

    if choices['is_purging'].get() == 0:
        type_purge = "no_purge"
    else:
        type_purge = "periodic_purge"

    if choices['is_precise'].get() == 0:
        max_step = np.inf
    else:
        max_step = 0.1

    if choices['is_synthetic'].get() == 0:
        type_display = "shattered_display"
    else:
        type_display = "synthetic_display"

    if choices['is_dynamic'].get() == 0:
        type_plot = "final_plot"
    else:
        type_plot = "dynamic_plot"

    return (Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, Aact, Hgdl, Hcl, Hmem, Hgc, Wgc, Lgc, epsilon_gdl,
            epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, t_step,
            i_step, i_max_pola, delta_pola, i_EIS, t_purge, delta_t_purge, type_fuel_cell, type_auxiliary, type_purge,
            max_step, type_display, type_plot)


def value_control(choices):
    """
    This function checks the integrity of the values entered by the user and 
    returns an empty tuple if they are not valid.
    """

    # The values entered by the user are checked for compliance
    if choices['Tfc_i'].get() < 0:
        messagebox.showerror(title='Temperatures', message=
        'Negative temperatures do not exist in the Kelvin scale.')
        choices.clear()
        return
    if choices['Pa_des_i'].get() < 0 or choices['Pc_des_i'].get() < 0 or \
            choices['Pc_des_i'].get() > 5.0 or choices['Pc_des_i'].get() > 5.0:
        messagebox.showerror(title='Desired pressures', message=
        'Desired pressure should be positive and bellow 5.0 bars.')
        choices.clear()
        return
    if choices['Sa_i'].get() < 1 or choices['Sa_i'].get() > 5 or \
            choices['Sc_i'].get() < 1 or choices['Sc_i'].get() > 5:
        messagebox.showerror(title='Stoichiometric ratios', message=
        'The stoichiometric ratios Sa and Sc should be between 1 and 5.')
        choices.clear()
        return
    if choices['Phi_a_des_i'].get() < 0 or choices['Phi_a_des_i'].get() > 1 or \
            choices['Phi_c_des_i'].get() < 0 or choices['Phi_c_des_i'].get() > 1:
        messagebox.showerror(title='Desired humidity', message=
        'The desired humidities should be between 0 and 1.')
        choices.clear()
        return
    if choices['Aact_i'].get() < 0:
        messagebox.showerror(title='Active area', message=
        'Negative active area is impossible.')
        choices.clear()
        return
    if choices['Hgdl_i'].get() < 1 or choices['Hgdl_i'].get() > 1000 or \
            choices['Hcl_i'].get() < 1 or choices['Hcl_i'].get() > 1000 or \
            choices['Hmem_i'].get() < 1 or choices['Hmem_i'].get() > 1000:
        messagebox.showerror(title='MEA thickness', message=
        'All MEA components generally have a thickness between 1µm and 1mm.')
        choices.clear()
        return
    if choices['Hgc_i'].get() < 10 or choices['Hgc_i'].get() > 10000 or \
            choices['Wgc_i'].get() < 10 or choices['Wgc_i'].get() > 10000 or \
            choices['Lgc_i'].get() < 0 or choices['Lgc_i'].get() > 100:
        messagebox.showerror(title='GC distances', message=
        'GC generally have a thickness and a width between 10µm and 10mm.\
        Also, GC length is generally between 0 and 100m')
        choices.clear()
        return
    if choices['epsilon_gdl_i'].get() < 0 or choices['epsilon_gdl_i'].get() > 1 or \
            choices['epsilon_mc_i'].get() < 0 or choices['epsilon_mc_i'].get() > 1:
        messagebox.showerror(title='Porosities', message=
        'All porosities should be between 0 and 1.')
        choices.clear()
        return
    if choices['tau_i'].get() < 1 or choices['tau_i'].get() > 4:
        messagebox.showerror(title='Pore structure coefficient', message=
        'The pore structure coefficient should be between 1 and 4.')
        choices.clear()
        return
    if choices['epsilon_c_i'].get() < 0 or choices['epsilon_c_i'].get() > 1:
        messagebox.showerror(title='Compression ratio', message=
        'The compression ratio should be between 0 and 1.')
        choices.clear()
        return
    if choices['e_i'].get() < 3 or choices['e_i'].get() > 5:
        messagebox.showerror(title='Capillary exponent', message=
        'The capillary exponent should be between 3 and 5 and being an integer.')
        choices.clear()
        return
    if choices['Re_i'].get() < 0.5 or choices['Re_i'].get() > 5:
        messagebox.showerror(title='Electron conduction resistance', message=
        'The electron conduction resistance is generally between 0.5 and 5 µΩ.m².')
        choices.clear()
        return
    if choices['i0_c_ref_i'].get() < 0.001 or choices['i0_c_ref_i'].get() > 500:
        messagebox.showerror(title='Referenced exchange current density', message=
        'The referenced exchange current density is generally between 0.001 and 500 A.m-2.')
        choices.clear()
        return
    if choices['kappa_co_i'].get() < 0.01 or choices['kappa_co_i'].get() > 100:
        messagebox.showerror(title='Crossover correction coefficient', message=
        'The crossover correction coefficient is generally between 0.01 and 100 mol.m-1.s-1.Pa-1.')
        choices.clear()
        return
    if choices['kappa_c_i'].get() < 0 or choices['kappa_c_i'].get() > 100:
        messagebox.showerror(title='Overpotential correction exponent', message=
        'The overpotential correction exponent is generally between 0 and 100.')
        choices.clear()
        return
    if choices['a_slim_i'].get() < 0 or choices['a_slim_i'].get() > 1:
        messagebox.showerror(title='Slop of slim function', message=
        'The slop of slim function is generally between 0 and 1.')
        choices.clear()
        return
    if choices['b_slim_i'].get() < 0 or choices['b_slim_i'].get() > 1:
        messagebox.showerror(title='Intercept of slim function', message=
        'The intercept of slim function is generally between 0 and 1.')
        choices.clear()
        return
    if choices['a_switch_i'].get() < 0 or choices['a_switch_i'].get() > 1:
        messagebox.showerror(title='Slop of switch function', message=
        'The slop of switch function is generally between 0 and 1.')
        choices.clear()
        return
    if choices['C_dl_i'].get() < 5 or choices['C_dl_i'].get() > 100:
        messagebox.showerror(title='Double layer capacitance', message=
        'I have not settled yet a range for C_dl.')
        choices.clear()
        return
    if choices['t0_step_i'].get() < 0 or choices['tf_step_i'].get() < 0 or choices['delta_t_load_step_i'].get() < 0 or \
            choices['delta_t_dyn_step_i'].get() < 0 or choices['delta_t_load_pola_i'].get() < 0 or \
            choices['delta_t_break_pola_i'].get() < 0 or choices['delta_t_ini_pola_i'].get() < 0 or \
            choices['t0_step_i'].get() > choices['tf_step_i'].get() or \
            choices['delta_t_load_step_i'].get() > (choices['tf_step_i'].get() - choices['t0_step_i'].get()):
        messagebox.showerror(title='Times', message=
        'The times should be positive, t0_step < tf_step and delta_t_load_step < (tf_step - t0_step).')
        choices.clear()
        return
    if choices['i_ini_step_i'].get() < 0 or choices['i_final_step_i'].get() < 0 or \
            choices['i_max_pola_i'].get() < 0 or choices['delta_i_pola_i'].get() < 0 or \
            choices['i_EIS_i'].get() < 0 or \
            choices['delta_i_pola_i'].get() > choices['i_max_pola_i'].get() or \
            choices['i_ini_step_i'].get() > choices['i_final_step_i'].get():
        messagebox.showerror(title='Current densities', message=
        'The current densities should be positive, delta_i_pola_i < i_max_pola_i and i_ini_step_i < i_final_step_i.')
        choices.clear()
        return
    if choices['t_purge_i'].get() < 0 and choices['delta_purge_i'].get() < 0:
        messagebox.showerror(title='Purge times', message=
        'Negative times does not characterise purges.')
        choices.clear()
        return
