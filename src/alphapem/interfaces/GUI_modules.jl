# -*- coding: utf-8 -*-

"""This module contains some of the required functions for the GUI.jl file.
"""

module GUIModules

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using PyCall
import Base.time as _time

# Importing Python modules for tkinter and message boxes
const tk = PyCall.pyimport("tkinter")
const messagebox = PyCall.pyimport("tkinter.messagebox")
const ttk = PyCall.pyimport("tkinter.ttk")
const plt = PyCall.pyimport("matplotlib.pyplot")
const np = PyCall.pyimport("numpy")

# Importing Julia modules
# using AlphaPEM.Config.Parameters: calculate_current_density_parameters, calculate_computing_parameters
# using AlphaPEM.Config.CurrentDensities: EIS_parameters
# using AlphaPEM.Config.ParametersSpecific: stored_operating_inputs, stored_physical_parameters
# using AlphaPEM.Application.RunSimulationModules: figures_preparation

# _____________________________________________________GUI modules_____________________________________________________

"""
    changeValue(...)

This function is called when the user selects a specific option from a dropdown menu for the type of fuel cell.
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
function changeValue(operating_conditions_frame, accessible_parameters_frame, undetermined_parameters_frame,
                     current_density_parameters_frame, computing_parameters_frame, choice_operating_conditions,
                     choice_accessible_parameters, choice_undetermined_parameters, choice_current_density_parameters,
                     choice_computing_parameters, choices_buttons)

    if choices_buttons["type_fuel_cell"]["value"].get() != "Enter your specifications"
        # Recovers the new settings
        recover_for_display_operating_inputs_and_physical_parameters(choice_operating_conditions,
                                                                     choice_accessible_parameters,
                                                                     choice_undetermined_parameters,
                                                                     choice_current_density_parameters,
                                                                     choice_computing_parameters, choices_buttons)
        # Display the labels for ...
        #       operating conditions
        for (k, v) in choice_operating_conditions
            ttk.Label(operating_conditions_frame, width=7, anchor="w", textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       accessible physical parameters
        for (k, v) in choice_accessible_parameters
            ttk.Label(accessible_parameters_frame, width=7, anchor="w", textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       undetermined physical parameters
        for (k, v) in choice_undetermined_parameters
            ttk.Label(undetermined_parameters_frame, width=7, anchor="w", textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       current density parameters
        for (k, v) in choice_current_density_parameters
            ttk.Label(current_density_parameters_frame, width=7, anchor="w", textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       computing parameters
        for (k, v) in choice_computing_parameters
            ttk.Label(computing_parameters_frame, width=7, anchor="w", textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end

    else  # choices_buttons["type_fuel_cell"]["value"].get() == "Enter your specifications"
        # Saves and displays the user entries for ...
        #       operating conditions
        for (k, v) in choice_operating_conditions
            ttk.Entry(operating_conditions_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       accessible physical parameters
        for (k, v) in choice_accessible_parameters
            ttk.Entry(accessible_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       undetermined physical parameters
        for (k, v) in choice_undetermined_parameters
            ttk.Entry(undetermined_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       current density parameters
        for (k, v) in choice_current_density_parameters
            ttk.Entry(current_density_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
        #       computing parameters
        for (k, v) in choice_computing_parameters
            ttk.Entry(computing_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
        end
    end
end

"""
    display_parameter_labels(...)

This function displays labels on the GUI, representing operating conditions and physical parameters, without
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
function display_parameter_labels(operating_conditions_frame, accessible_parameters_frame, undetermined_parameters_frame,
                                  current_density_parameters_frame, computing_parameters_frame, choice_operating_conditions,
                                  choice_accessible_parameters, choice_undetermined_parameters,
                                  choice_current_density_parameters, choice_computing_parameters)

    # Display the titles
    ttk.Label(operating_conditions_frame, text="Operating conditions", font=("cmr10", 12, "bold")).grid(row=1, column=0, columnspan=6, ipady=15)
    ttk.Label(accessible_parameters_frame, text="Accessible physical parameters", font=("cmr10", 12, "bold")).grid(row=0, column=0, columnspan=6, ipady=15)

    # Display the labels for ...
    #       operating conditions
    for (k, v) in choice_operating_conditions
        ttk.Label(operating_conditions_frame, text=k, font=("cmr10", 10)).grid(row=v["label_row"], column=v["label_column"] - 1, sticky="w")
    end
    #       accessible physical parameters
    for (k, v) in choice_accessible_parameters
        ttk.Label(accessible_parameters_frame, text=k, font=("cmr10", 10)).grid(row=v["label_row"], column=v["label_column"] - 1, sticky="w")
    end
    #       undetermined physical parameters
    for (k, v) in choice_undetermined_parameters
        ttk.Label(undetermined_parameters_frame, text=k, font=("cmr10", 10)).grid(row=v["label_row"], column=v["label_column"] - 1, sticky="w")
    end
    #       current density parameters
    ttk.Label(current_density_parameters_frame, text="Step current parameters", font=("cmr10", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w")
    ttk.Label(current_density_parameters_frame, text="Polarization current parameters", font=("cmr10", 10, "bold")).grid(row=2, column=0, columnspan=2, sticky="w")
    ttk.Label(current_density_parameters_frame, text="EIS current parameters", font=("cmr10", 10, "bold")).grid(row=5, column=0, columnspan=2, sticky="w")
    for (k, v) in choice_current_density_parameters
        ttk.Label(current_density_parameters_frame, text=k, font=("cmr10", 10)).grid(row=v["label_row"], column=v["label_column"] - 1, sticky="w")
    end
    #       computing parameters
    for (k, v) in choice_computing_parameters
        ttk.Label(computing_parameters_frame, text=k, font=("cmr10", 10)).grid(row=v["label_row"], column=v["label_column"] - 1, sticky="w")
    end
end

"""
    display_parameters_value(...)

This function displays entry widgets on the GUI. There, the user can enter values for operating conditions and
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
function display_parameters_value(operating_conditions_frame, accessible_parameters_frame, undetermined_parameters_frame,
                                  current_density_parameters_frame, computing_parameters_frame, choice_operating_conditions,
                                  choice_accessible_parameters, choice_undetermined_parameters,
                                  choice_current_density_parameters, choice_computing_parameters)

    # Display the value for ...
    #       operating conditions
    for (k, v) in choice_operating_conditions
        ttk.Entry(operating_conditions_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
    end
    #       accessible physical parameters
    for (k, v) in choice_accessible_parameters
        ttk.Entry(accessible_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
    end
    #       undetermined physical parameters
    for (k, v) in choice_undetermined_parameters
        ttk.Entry(undetermined_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
    end
    #       current density parameters
    for (k, v) in choice_current_density_parameters
        ttk.Entry(current_density_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
    end
    #       computing parameters
    for (k, v) in choice_computing_parameters
        ttk.Entry(computing_parameters_frame, width=7, textvariable=v["value"]).grid(row=v["label_row"], column=v["label_column"], padx=5)
    end
end

"""
    display_radiobuttons(...)

This function displays radiobuttons on the GUI, allowing the user to make choices for control, results display,
plot style, etc.

Parameters
----------
model_possibilities_frame : ttk.Frame
    The frame where the graphical elements for the model possibilities and the choice of current density are placed.
choices_buttons : dict
    A dictionary containing the button information.
"""
function display_radiobuttons(model_possibilities_frame, choices_buttons)

    ttk.Label(model_possibilities_frame, text="Model configuration", font=("cmr10", 12, "bold")).grid(row=0, column=0, columnspan=6, ipady=15)

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text="Auxiliaries: ", font=("cmr10", 12)).grid(row=choices_buttons["type_auxiliary"]["label_row"], column=0, columnspan=1, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="No auxiliaries", value=0,
                    variable=choices_buttons["type_auxiliary"]["value"]).grid(row=choices_buttons["type_auxiliary"]["label_row"], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Forced-convective cathode\nwith anodic recirculation", value=1,
                    variable=choices_buttons["type_auxiliary"]["value"]).grid(row=choices_buttons["type_auxiliary"]["label_row"], column=3, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Forced-convective cathode\nwith flow-through anode", value=2,
                    variable=choices_buttons["type_auxiliary"]["value"]).grid(row=choices_buttons["type_auxiliary"]["label_row"], column=4, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text="Voltage zone: ", font=("cmr10", 12)).grid(row=choices_buttons["voltage_zone"]["label_row"], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Full", value=0,
                    variable=choices_buttons["voltage_zone"]["value"]).grid(row=choices_buttons["voltage_zone"]["label_row"], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Before voltage drop", value=1,
                    variable=choices_buttons["voltage_zone"]["value"]).grid(row=choices_buttons["voltage_zone"]["label_row"], column=3, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text="Purge: ", font=("cmr10", 12)).grid(row=choices_buttons["type_purge"]["label_row"], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="No purge", value=0,
                    variable=choices_buttons["type_purge"]["value"]).grid(row=choices_buttons["type_purge"]["label_row"], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Periodic", value=1,
                    variable=choices_buttons["type_purge"]["value"]).grid(row=choices_buttons["type_purge"]["label_row"], column=3, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Constant", value=2,
                    variable=choices_buttons["type_purge"]["value"]).grid(row=choices_buttons["type_purge"]["label_row"], column=4, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text="Display: ", font=("cmr10", 12)).grid(row=choices_buttons["type_display"]["label_row"], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="No display", value=0,
                    variable=choices_buttons["type_display"]["value"]).grid(row=choices_buttons["type_display"]["label_row"], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Synthetic", value=1,
                    variable=choices_buttons["type_display"]["value"]).grid(row=choices_buttons["type_display"]["label_row"], column=3, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Multiple", value=2,
                    variable=choices_buttons["type_display"]["value"]).grid(row=choices_buttons["type_display"]["label_row"], column=4, sticky="w")

    # Ask the user to choose an option and save it
    ttk.Label(model_possibilities_frame, text="Plot: ", font=("cmr10", 12)).grid(row=choices_buttons["type_plot"]["label_row"], column=0, columnspan=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Fixed", value=0,
                    variable=choices_buttons["type_plot"]["value"]).grid(row=choices_buttons["type_plot"]["label_row"], column=2, sticky="w")
    ttk.Radiobutton(model_possibilities_frame, text="Dynamic", value=1,
                    variable=choices_buttons["type_plot"]["value"]).grid(row=choices_buttons["type_plot"]["label_row"], column=3, sticky="w")
end

"""
    recover_for_display_operating_inputs_and_physical_parameters(...)

This function retrieves parameter values for predefined stacks (e.g., "EH-31 1.5 bar (2021)", "Biao Xie 1.0 bar
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
function recover_for_display_operating_inputs_and_physical_parameters(choice_operating_conditions,
                                                                      choice_accessible_parameters,
                                                                      choice_undetermined_parameters,
                                                                      choice_current_density_parameters,
                                                                      choice_computing_parameters, choice_buttons)

    # type_fuel_cell recovery
    if choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack (2022)"
        type_fuel_cell = "ZSW-GenStack"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack: Pa=1.61 bar, Pc=1.41 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_1.61_Pc_1.41"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack: Pa=2.01 bar, Pc=1.81 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_2.01_Pc_1.81"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack: Pa=2.4 bar, Pc=2.2 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_2.4_Pc_2.2"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack: Pa=2.8 bar, Pc=2.6 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_2.8_Pc_2.6"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack: T=62°C (2022)"
        type_fuel_cell = "ZSW-GenStack_T_62"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack: T=76°C (2022)"
        type_fuel_cell = "ZSW-GenStack_T_76"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "ZSW-GenStack: T=84°C (2022)"
        type_fuel_cell = "ZSW-GenStack_T_84"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "EH-31: P=1.5 bar (2021)"
        type_fuel_cell = "EH-31_1.5"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "EH-31: P=2.0 bar (2021)"
        type_fuel_cell = "EH-31_2.0"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "EH-31: P=2.25 bar (2021)"
        type_fuel_cell = "EH-31_2.25"
    elseif choice_buttons["type_fuel_cell"]["value"].get() == "EH-31: P=2.5 bar (2021)"
        type_fuel_cell = "EH-31_2.5"
    end

    if choice_buttons["voltage_zone"]["value"].get() == 0
        voltage_zone = "full"
    else
        voltage_zone = "before_voltage_drop"
    end

    # ... Rest of the function would be uncommented when implementation is complete
    # (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
    # current_density) = calculate_current_density_parameters()
    #
    # T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, i_max_pola = stored_operating_inputs(type_fuel_cell, voltage_zone)
    #
    # (Hacl, Hccl, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc, Lgc, nb_channel_in_gc,
    #  Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem, Aact, nb_cell, e, K_l_ads, K_O2_ad_Pt, Re, i0_c_ref, kappa_co,
    #  kappa_c, C_scl) = stored_physical_parameters(type_fuel_cell)
    #
    # nb_gc, nb_gdl, nb_mpl, t_purge, rtol, atol = calculate_computing_parameters(step_current_parameters)
    #
    # # operating conditions recovery
    # choice_operating_conditions["Temperature - Tfc (°C)"]["value"].set(round(T_des - 273.15, 4))  # °C
    # ...
end

"""
    recover_for_use_operating_inputs_and_physical_parameters(...)

This function retrieves and converts the parameter values from the GUI into standard units for further
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
function recover_for_use_operating_inputs_and_physical_parameters(choice_operating_conditions, choice_accessible_parameters,
                                                                  choice_undetermined_parameters,
                                                                  choice_current_density_parameters,
                                                                  choice_computing_parameters, choice_buttons)
    # operating conditions
    T_des = choice_operating_conditions["Temperature - Tfc (°C)"]["value"].get() + 273.15  # K
    Pa_des = choice_operating_conditions["Anode pressure - Pa (bar)"]["value"].get() * 1e5  # Pa
    Pc_des = choice_operating_conditions["Cathode pressure - Pc (bar)"]["value"].get() * 1e5  # Pa
    Sa = choice_operating_conditions["Anode stoichiometry - Sa"]["value"].get()
    Sc = choice_operating_conditions["Cathode stoichiometry - Sc"]["value"].get()
    Phi_a_des = choice_operating_conditions["Anode humidity - Φa"]["value"].get()
    Phi_c_des = choice_operating_conditions["Cathode humidity - Φc"]["value"].get()
    y_H2_in = choice_operating_conditions["Anode inlet H2 ratio - y_H2_in\n(flow-through anode only)"]["value"].get()
    # accessible physical parameters
    Aact = choice_accessible_parameters["Active area - Aact (cm²)"]["value"].get() * 1e-4  # m²
    nb_cell = Int(choice_accessible_parameters["Number of cells - nb_cell"]["value"].get())
    Hagc = choice_accessible_parameters["Anode gas channel\nthickness - Hagc (µm)"]["value"].get() * 1e-6  # m
    Hcgc = choice_accessible_parameters["Cathode gas channel\nthickness - Hcgc (µm)"]["value"].get() * 1e-6  # m
    Wagc = choice_accessible_parameters["Anode gas channel\nwidth - Wagc (µm)"]["value"].get() * 1e-6  # m
    Wcgc = choice_accessible_parameters["Cathode gas channel\nwidth - Wcgc (µm)"]["value"].get() * 1e-6  # m
    Lgc = choice_accessible_parameters["Gas channel length - Lgc (mm)"]["value"].get() * 1e-3  # m
    nb_channel_in_gc = choice_accessible_parameters["Number of channels inside the\ngas channel - nb_channel_in_gc"]["value"].get()
    Lm = choice_accessible_parameters["Manifold length - Lm (mm)"]["value"].get() * 1e-3  # m
    Ldist = choice_accessible_parameters["Distributor length - Ldist (mm)"]["value"].get() * 1e-3  # m
    A_T_a = choice_accessible_parameters["Exhaust anode manifold throttle\narea - A_T_a (cm²)"]["value"].get() * 1e-4  # m²
    A_T_c = choice_accessible_parameters["Exhaust cathode manifold throttle\narea - A_T_c (cm²)"]["value"].get() * 1e-4  # m²
    Vasm = choice_accessible_parameters["Supply anode manifold\nvolume - Vasm (cm³)"]["value"].get() * 1e-6  # m³
    Vcsm = choice_accessible_parameters["Supply cathode manifold\nvolume - Vcsm (cm³)"]["value"].get() * 1e-6  # m³
    Vaem = choice_accessible_parameters["Exhaust anode manifold\nvolume - Vaem (cm³)"]["value"].get() * 1e-6  # m³
    Vcem = choice_accessible_parameters["Exhaust cathode manifold\nvolume - Vcem (cm³)"]["value"].get() * 1e-6  # m³
    V_endplate_a = choice_accessible_parameters["Anode endplate volume - V_endplate_c (cm³)"]["value"].get() * 1e-6  # m³
    V_endplate_c = choice_accessible_parameters["Cathode endplate volume - V_endplate_c (cm³)"]["value"].get() * 1e-6  # m³
    # undetermined physical parameters
    Hgdl = choice_undetermined_parameters["GDL thickness - Hgdl (µm)"]["value"].get() * 1e-6  # m
    Hmpl = choice_undetermined_parameters["MPL thickness - Hmpl (µm)"]["value"].get() * 1e-6  # m
    Hacl = choice_undetermined_parameters["ACL thickness - Hacl (µm)"]["value"].get() * 1e-6  # m
    Hccl = choice_undetermined_parameters["CCL thickness - Hccl (µm)"]["value"].get() * 1e-6  # m
    Hmem = choice_undetermined_parameters["Membrane thickness - Hmem (µm)"]["value"].get() * 1e-6  # m
    epsilon_gdl = choice_undetermined_parameters["GDL porosity - ε_gdl"]["value"].get()
    epsilon_mpl = choice_undetermined_parameters["MPL porosity - ε_mpl"]["value"].get()
    epsilon_c = choice_undetermined_parameters["Compression ratio - ε_c"]["value"].get()
    e = choice_undetermined_parameters["Capillary exponent - e"]["value"].get()
    K_l_ads = choice_undetermined_parameters["Ratio of liq. and vap. water sorption\nrates in the mem. - K_l_ads"]["value"].get()
    K_O2_ad_Pt = choice_undetermined_parameters["Interfacial res. coef. of\nO2 ads. on Pt - K_O2_ad_Pt"]["value"].get()
    Re = choice_undetermined_parameters["Electron conduction\nresistance - Re (Ω.mm²)"]["value"].get() * 1e-6  # Ω.m²
    i0_c_ref = choice_undetermined_parameters["Reference exchange current\ndensity - i0_c_ref (A/m²)"]["value"].get()  # A.m-2
    kappa_co = choice_undetermined_parameters["Crossover correction coefficient\n- κ_co (mol/(m.s.Pa))"]["value"].get()  # mol.m-1.s-1.Pa-1
    kappa_c = choice_undetermined_parameters["Overpotential correction\nexponent - κ_c"]["value"].get()
    C_scl = choice_undetermined_parameters["Volumetric space-charge layer\ncapacitance - C_scl (F/cm³)"]["value"].get() * 1e6  # F.m-3

    # current density parameters
    delta_t_ini_step = choice_current_density_parameters["Stabilisation time\n- Δt_ini_step (min)"]["value"].get() * 60  # s
    delta_t_load_step = choice_current_density_parameters["Loading time\n- Δt_load_step (s)"]["value"].get()  # s
    delta_t_break_step = choice_current_density_parameters["Breaking time\n- Δt_break_step (min)"]["value"].get() * 60  # s
    i_step = choice_current_density_parameters["Current density step\n- i_step (A/cm²)"]["value"].get() * 1e4  # A.m-2
    delta_t_dyn_step = choice_computing_parameters["Time for dynamic\ndisplay - Δt_dyn_step (s)"]["value"].get()  # s
    step_current_parameters = Dict("delta_t_ini_step" => delta_t_ini_step, "delta_t_load_step" => delta_t_load_step,
                                   "delta_t_break_step" => delta_t_break_step, "i_step" => i_step,
                                   "delta_t_dyn_step" => delta_t_dyn_step)

    delta_t_ini_pola = choice_current_density_parameters["Stabilisation time\n- Δt_ini_pola (min)"]["value"].get() * 60  # s
    delta_t_load_pola = choice_current_density_parameters["Loading time\n- Δt_load_pola (s)"]["value"].get()  # s
    delta_t_break_pola = choice_current_density_parameters["Breaking time\n- Δt_break_pola (min)"]["value"].get() * 60  # s
    delta_i_pola = choice_current_density_parameters["Current density step\n- Δi_pola (A/cm²)"]["value"].get() * 1e4  # A.m-2
    i_max_pola = choice_current_density_parameters["Maximum current density\n- i_max_pola (A/cm²)"]["value"].get() * 1e4  # A.m-2
    pola_current_parameters = Dict("delta_t_ini_pola" => delta_t_ini_pola, "delta_t_load_pola" => delta_t_load_pola,
                                   "delta_t_break_pola" => delta_t_break_pola, "delta_i_pola" => delta_i_pola,
                                   "i_max_pola" => i_max_pola)

    pola_current_for_cali_parameters = nothing  # Calibration is not implemented in the GUI.
    i_EIS = choice_current_density_parameters["Static current\n- i_EIS (A/cm²)"]["value"].get() * 1e4  # (A.m-2)
    ratio_EIS = choice_current_density_parameters["Current ratio\n- ratio_EIS (%)"]["value"].get() / 100
    f_EIS = (choice_current_density_parameters["Power of the\ninitial frequency\n- f_power_min_EIS"]["value"].get(),
             choice_current_density_parameters["Power of the\nfinal frequency\n- f_power_max_EIS"]["value"].get(),
             choice_current_density_parameters["Number of frequencies\ntested - nb_f_EIS"]["value"].get(),
             choice_current_density_parameters["Number of points\ncalculated - nb_points_EIS"]["value"].get())
    # t_EIS = EIS_parameters(f_EIS)  # Time parameters for the EIS_current density function.
    t_EIS = nothing  # Placeholder - would call EIS_parameters when implementation is complete

    # computing parameters
    t_purge = pycall(choice_computing_parameters["Purge time - t_purge (s)"]["value"].get, PyAny)  # s
    delta_t_purge = pycall(choice_computing_parameters["Time between two purges\n- Δt_purge (s)"]["value"].get, PyAny)  # s
    nb_gc = pycall(choice_computing_parameters["Number of GC nodes - nb_gc"]["value"].get, PyAny)
    nb_gdl = pycall(choice_computing_parameters["Number of GDL nodes - nb_gdl"]["value"].get, PyAny)
    nb_mpl = pycall(choice_computing_parameters["Number of MPL nodes - nb_mpl"]["value"].get, PyAny)
    rtol = pycall(choice_computing_parameters["Solver relative tolerance - rtol"]["value"].get, PyAny)
    atol = pycall(choice_computing_parameters["Solver absolute tolerance - atol"]["value"].get, PyAny)

    fuel_cell_value = pycall(choice_buttons["type_fuel_cell"]["value"].get, PyAny)
    if fuel_cell_value == "ZSW-GenStack (2022)"
        type_fuel_cell = "ZSW-GenStack"
    elseif fuel_cell_value == "ZSW-GenStack: Pa=1.61 bar, Pc=1.41 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_1.61_Pc_1.41"
    elseif fuel_cell_value == "ZSW-GenStack: Pa=2.01 bar, Pc=1.81 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_2.01_Pc_1.81"
    elseif fuel_cell_value == "ZSW-GenStack: Pa=2.4 bar, Pc=2.2 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_2.4_Pc_2.2"
    elseif fuel_cell_value == "ZSW-GenStack: Pa=2.8 bar, Pc=2.6 bar (2022)"
        type_fuel_cell = "ZSW-GenStack_Pa_2.8_Pc_2.6"
    elseif fuel_cell_value == "ZSW-GenStack: T=62°C (2022)"
        type_fuel_cell = "ZSW-GenStack_T_62"
    elseif fuel_cell_value == "ZSW-GenStack: T=76°C (2022)"
        type_fuel_cell = "ZSW-GenStack_T_76"
    elseif fuel_cell_value == "ZSW-GenStack: T=84°C (2022)"
        type_fuel_cell = "ZSW-GenStack_T_84"
    elseif fuel_cell_value == "EH-31: P=1.5 bar (2021)"
        type_fuel_cell = "EH-31_1.5"
    elseif fuel_cell_value == "EH-31: P=2.0 bar (2021)"
        type_fuel_cell = "EH-31_2.0"
    elseif fuel_cell_value == "EH-31: P=2.25 bar (2021)"
        type_fuel_cell = "EH-31_2.25"
    elseif fuel_cell_value == "EH-31: P=2.5 bar (2021)"
        type_fuel_cell = "EH-31_2.5"
    elseif fuel_cell_value == "Enter your specifications"
        type_fuel_cell = "manual_setup"
    end

    auxiliary_value = Int(pycall(choice_buttons["type_auxiliary"]["value"].get, PyAny))
    if auxiliary_value == 0
        type_auxiliary = "no_auxiliary"
    elseif auxiliary_value == 1
        type_auxiliary = "forced-convective_cathode_with_anodic_recirculation"
    else
        type_auxiliary = "forced-convective_cathode_with_flow-through_anode"
    end

    voltage_zone_value = Int(pycall(choice_buttons["voltage_zone"]["value"].get, PyAny))
    if voltage_zone_value == 0
        voltage_zone = "full"
    else
        voltage_zone = "before_voltage_drop"
    end

    purge_value = Int(pycall(choice_buttons["type_purge"]["value"].get, PyAny))
    if purge_value == 0
        type_purge = "no_purge"
    elseif purge_value == 1
        type_purge = "periodic_purge"
    else
        type_purge = "constant_purge"
    end

    display_value = Int(pycall(choice_buttons["type_display"]["value"].get, PyAny))
    if display_value == 0
        type_display = "no_display"
    elseif display_value == 1
        type_display = "synthetic"
    else
        type_display = "multiple"
    end

    if choice_buttons["type_plot"]["value"].get() == 0
        type_plot = "fixed"
    else
        type_plot = "dynamic"
    end

    return (T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, Aact, nb_cell, Hgdl, Hmpl, Hacl, Hccl, Hmem,
            Hagc, Hcgc, Wagc, Wcgc, Lgc, nb_channel_in_gc, Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem,
            V_endplate_a, V_endplate_c, epsilon_gdl, epsilon_mpl, epsilon_c, e, K_l_ads, K_O2_ad_Pt, Re, i0_c_ref,
            kappa_co, kappa_c, C_scl, step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters,
            i_EIS, ratio_EIS, f_EIS, t_EIS, t_purge, delta_t_purge, nb_gc, nb_gdl, nb_mpl, rtol, atol, type_fuel_cell,
            voltage_zone, type_auxiliary, type_purge, type_display, type_plot)
end

"""
    value_control(...)

This function checks the integrity of the values entered by the user and returns an empty tuple if they are not
valid.

Parameters
----------
choice_operating_conditions : dict
    A dictionary containing the operating condition information.
choice_accessible_parameters : dict
    A dictionary containing the accessible parameter information.
choice_undetermined_parameters : dict
    A dictionary containing the undetermined parameter information.
choice_current_density_parameters : dict
    A dictionary containing the current density parameter information.
choice_computing_parameters : dict
    A dictionary containing the computing parameter information.
choice_buttons : dict
    A dictionary containing the button information.
current_button : int
    A dictionary representing the clicked button.
"""
function value_control(choice_operating_conditions, choice_accessible_parameters, choice_undetermined_parameters,
                       choice_current_density_parameters, choice_computing_parameters, choice_buttons, current_button)

    # The values entered by the user are checked for compliance
    if choice_operating_conditions["Temperature - Tfc (°C)"]["value"].get() < 0
        messagebox.showerror(title="Temperatures", message="Negative temperatures do not exist in the Kelvin scale.")
        return
    end

    if choice_operating_conditions["Anode pressure - Pa (bar)"]["value"].get() < 0 ||
            choice_operating_conditions["Cathode pressure - Pc (bar)"]["value"].get() < 0 ||
            choice_operating_conditions["Cathode pressure - Pc (bar)"]["value"].get() > 5.0 ||
            choice_operating_conditions["Cathode pressure - Pc (bar)"]["value"].get() > 5.0
        messagebox.showerror(title="Desired pressures", message="Desired pressure should be positive and bellow 5.0 bars.")
        return
    end

    if choice_operating_conditions["Anode stoichiometry - Sa"]["value"].get() < 1 ||
            choice_operating_conditions["Anode stoichiometry - Sa"]["value"].get() > 5 ||
            choice_operating_conditions["Cathode stoichiometry - Sc"]["value"].get() < 1 ||
            choice_operating_conditions["Cathode stoichiometry - Sc"]["value"].get() > 5
        messagebox.showerror(title="Stoichiometric ratios", message="The stoichiometric ratios Sa and Sc should be between 1 and 5.")
        return
    end

    if choice_operating_conditions["Anode humidity - Φa"]["value"].get() < 0 ||
            choice_operating_conditions["Anode humidity - Φa"]["value"].get() > 1 ||
            choice_operating_conditions["Cathode humidity - Φc"]["value"].get() < 0 ||
            choice_operating_conditions["Cathode humidity - Φc"]["value"].get() > 1
        messagebox.showerror(title="Desired humidity", message="The desired humidities should be between 0 and 1.")
        return
    end

    if choice_operating_conditions["Anode inlet H2 ratio - y_H2_in\n(flow-through anode only)"]["value"].get() < 0 ||
            choice_operating_conditions["Anode inlet H2 ratio - y_H2_in\n(flow-through anode only)"]["value"].get() > 1
        messagebox.showerror(title="Anode inlet H2 ratio", message="The anode inlet H2 ratio should be between 0 and 1.")
        return
    end

    if choice_accessible_parameters["Active area - Aact (cm²)"]["value"].get() < 0
        messagebox.showerror(title="Active area", message="Negative active area is impossible.")
        return
    end

    if choice_accessible_parameters["Number of cells - nb_cell"]["value"].get() < 0
        messagebox.showerror(title="Number of cells", message="Negative number of cells is impossible.")
        return
    end

    if choice_accessible_parameters["Exhaust anode manifold throttle\narea - A_T_a (cm²)"]["value"].get() < 0 ||
            choice_accessible_parameters["Exhaust cathode manifold throttle\narea - A_T_c (cm²)"]["value"].get() < 0 ||
            choice_accessible_parameters["Supply anode manifold\nvolume - Vasm (cm³)"]["value"].get() < 0 ||
            choice_accessible_parameters["Supply cathode manifold\nvolume - Vcsm (cm³)"]["value"].get() < 0 ||
            choice_accessible_parameters["Exhaust anode manifold\nvolume - Vaem (cm³)"]["value"].get() < 0 ||
            choice_accessible_parameters["Exhaust cathode manifold\nvolume - Vcem (cm³)"]["value"].get() < 0 ||
            choice_accessible_parameters["Anode endplate volume - V_endplate_c (cm³)"]["value"].get() < 0 ||
            choice_accessible_parameters["Cathode endplate volume - V_endplate_c (cm³)"]["value"].get() < 0
        messagebox.showerror(title="Manifold parameters", message="Negative volumes or area are impossible.")
        return
    end

    if choice_accessible_parameters["Anode gas channel\nthickness - Hagc (µm)"]["value"].get() < 10 ||
            choice_accessible_parameters["Anode gas channel\nthickness - Hagc (µm)"]["value"].get() > 10000 ||
            choice_accessible_parameters["Cathode gas channel\nthickness - Hcgc (µm)"]["value"].get() < 10 ||
            choice_accessible_parameters["Cathode gas channel\nthickness - Hcgc (µm)"]["value"].get() > 10000 ||
            choice_accessible_parameters["Anode gas channel\nwidth - Wagc (µm)"]["value"].get() < 10 ||
            choice_accessible_parameters["Anode gas channel\nwidth - Wagc (µm)"]["value"].get() > 10000 ||
            choice_accessible_parameters["Cathode gas channel\nwidth - Wcgc (µm)"]["value"].get() < 10 ||
            choice_accessible_parameters["Cathode gas channel\nwidth - Wcgc (µm)"]["value"].get() > 10000 ||
            choice_accessible_parameters["Gas channel length - Lgc (mm)"]["value"].get() < 0 ||
            choice_accessible_parameters["Gas channel length - Lgc (mm)"]["value"].get() > 1000 ||
            choice_accessible_parameters["Distributor length - Ldist (mm)"]["value"].get() < 0 ||
            choice_accessible_parameters["Manifold length - Lm (mm)"]["value"].get() < 0
        messagebox.showerror(title="GC distances", message="GC generally have a thickness and a width between 10µm and 10mm. Also, GC length is generally between 0 and 1000mm")
        return
    end

    if choice_accessible_parameters["Number of channels inside the\ngas channel - nb_channel_in_gc"]["value"].get() < 1 ||
            typeof(choice_accessible_parameters["Number of channels inside the\ngas channel - nb_channel_in_gc"]["value"].get()) != Int
        messagebox.showerror(title="nb_channel_in_gc", message="The nb_channel_in_gc value should be an integer bigger or equal to 1.")
        return
    end

    if choice_undetermined_parameters["GDL thickness - Hgdl (µm)"]["value"].get() < 1 ||
            choice_undetermined_parameters["GDL thickness - Hgdl (µm)"]["value"].get() > 1000 ||
            choice_undetermined_parameters["MPL thickness - Hmpl (µm)"]["value"].get() < 1 ||
            choice_undetermined_parameters["MPL thickness - Hmpl (µm)"]["value"].get() > 1000 ||
            choice_undetermined_parameters["ACL thickness - Hacl (µm)"]["value"].get() < 1 ||
            choice_undetermined_parameters["ACL thickness - Hacl (µm)"]["value"].get() > 1000 ||
            choice_undetermined_parameters["CCL thickness - Hccl (µm)"]["value"].get() < 1 ||
            choice_undetermined_parameters["CCL thickness - Hccl (µm)"]["value"].get() > 1000 ||
            choice_undetermined_parameters["Membrane thickness - Hmem (µm)"]["value"].get() < 1 ||
            choice_undetermined_parameters["Membrane thickness - Hmem (µm)"]["value"].get() > 1000
        messagebox.showerror(title="MEA thickness", message="All MEA components generally have a thickness between 1µm and 1mm.")
        return
    end

    if choice_undetermined_parameters["GDL porosity - ε_gdl"]["value"].get() < 0.4 ||
            choice_undetermined_parameters["GDL porosity - ε_gdl"]["value"].get() > 0.95
        messagebox.showerror(title="GDL porosity", message="GDL porosity should be between 0.4 and 0.95.")
        return
    end

    if choice_undetermined_parameters["MPL porosity - ε_mpl"]["value"].get() < 0.30 ||
            choice_undetermined_parameters["MPL porosity - ε_mpl"]["value"].get() > 0.70
        messagebox.showerror(title="MPL porosity", message="MPL porosity should be between 0.30 and 0.70.")
        return
    end

    if choice_undetermined_parameters["Compression ratio - ε_c"]["value"].get() < 0 ||
            choice_undetermined_parameters["Compression ratio - ε_c"]["value"].get() > 1
        messagebox.showerror(title="Compression ratio", message="The compression ratio should be between 0 and 1.")
        return
    end

    if choice_undetermined_parameters["Capillary exponent - e"]["value"].get() < 3 ||
            choice_undetermined_parameters["Capillary exponent - e"]["value"].get() > 5
        messagebox.showerror(title="Capillary exponent", message="The capillary exponent should be between 3 and 5 and being an integer.")
        return
    end

    if choice_undetermined_parameters["Ratio of liq. and vap. water sorption\nrates in the mem. - K_l_ads"]["value"].get() < 0 ||
            choice_undetermined_parameters["Ratio of liq. and vap. water sorption\nrates in the mem. - K_l_ads"]["value"].get() > 100
        messagebox.showerror(title="Ratio of liquid and vapor water sorption rates in the membrane",
                             message="The ratio of liquid and vapour water sorption rates in the membrane should be between 0 and 100.")
    end

    if choice_undetermined_parameters["Interfacial res. coef. of\nO2 ads. on Pt - K_O2_ad_Pt"]["value"].get() < 1 ||
            choice_undetermined_parameters["Interfacial res. coef. of\nO2 ads. on Pt - K_O2_ad_Pt"]["value"].get() > 10
        messagebox.showerror(title="Interfacial res. coef. of O2 ads. on Pt", message="The interfacial res. coef. of O2 ads. on Pt should be between 10 and 10")
        return
    end

    if choice_undetermined_parameters["Electron conduction\nresistance - Re (Ω.mm²)"]["value"].get() < 0
        messagebox.showerror(title="Electron conduction resistance", message="Re should be positive.")
        return
    end

    if choice_undetermined_parameters["Reference exchange current\ndensity - i0_c_ref (A/m²)"]["value"].get() < 0.001 ||
            choice_undetermined_parameters["Reference exchange current\ndensity - i0_c_ref (A/m²)"]["value"].get() > 100
        messagebox.showerror(title="Referenced exchange current densities", message="The referenced exchange current density is generally between 0.001 and 100 A.m-2.")
        return
    end

    if choice_undetermined_parameters["Crossover correction coefficient\n- κ_co (mol/(m.s.Pa))"]["value"].get() < 0.01 ||
            choice_undetermined_parameters["Crossover correction coefficient\n- κ_co (mol/(m.s.Pa))"]["value"].get() > 100
        messagebox.showerror(title="Crossover correction coefficient", message="The crossover correction coefficient is generally between 0.01 and 100 mol.m-1.s-1.Pa-1.")
        return
    end

    if choice_undetermined_parameters["Overpotential correction\nexponent - κ_c"]["value"].get() < 0 ||
            choice_undetermined_parameters["Overpotential correction\nexponent - κ_c"]["value"].get() > 100
        messagebox.showerror(title="Overpotential correction exponent", message="The overpotential correction exponent is generally between 0 and 100.")
        return
    end

    if choice_undetermined_parameters["Volumetric space-charge layer\ncapacitance - C_scl (F/cm³)"]["value"].get() < 5 ||
            choice_undetermined_parameters["Volumetric space-charge layer\ncapacitance - C_scl (F/cm³)"]["value"].get() > 100
        messagebox.showerror(title="Double layer capacitance", message="I have not settled yet a range for C_scl.")
        return
    end

    if choice_current_density_parameters["Stabilisation time\n- Δt_ini_step (min)"]["value"].get() < 0 ||
            choice_current_density_parameters["Loading time\n- Δt_load_step (s)"]["value"].get() < 0 ||
            choice_current_density_parameters["Breaking time\n- Δt_break_step (min)"]["value"].get() < 0 ||
            choice_computing_parameters["Time for dynamic\ndisplay - Δt_dyn_step (s)"]["value"].get() < 0 ||
            choice_current_density_parameters["Stabilisation time\n- Δt_ini_pola (min)"]["value"].get() < 0 ||
            choice_current_density_parameters["Loading time\n- Δt_load_pola (s)"]["value"].get() < 0 ||
            choice_current_density_parameters["Breaking time\n- Δt_break_pola (min)"]["value"].get() < 0
        messagebox.showerror(title="Times", message="The times should be positive, t0_step < tf_step and delta_t_load_step < (tf_step - t0_step).")
        return
    end

    if choice_current_density_parameters["Maximum current density\n- i_max_pola (A/cm²)"]["value"].get() < 0 ||
            choice_current_density_parameters["Current density step\n- Δi_pola (A/cm²)"]["value"].get() < 0 ||
            choice_current_density_parameters["Static current\n- i_EIS (A/cm²)"]["value"].get() < 0 ||
            choice_current_density_parameters["Current density step\n- Δi_pola (A/cm²)"]["value"].get() >
            choice_current_density_parameters["Maximum current density\n- i_max_pola (A/cm²)"]["value"].get()
        messagebox.showerror(title="Current densities", message="The current densities should be positive, delta_i_pola < i_max_pola and i_ini_step < i_final_step.")
        return
    end

    if choice_current_density_parameters["Current ratio\n- ratio_EIS (%)"]["value"].get() < 0 ||
            choice_current_density_parameters["Current ratio\n- ratio_EIS (%)"]["value"].get() > 20
        messagebox.showerror(title="Ratio EIS", message="Ratio EIS is a percentage of i_EIS and should be between 0 and 20 for plotting correct EIS.")
        return
    end

    if choice_current_density_parameters["Number of frequencies\ntested - nb_f_EIS"]["value"].get() < 0 ||
            choice_current_density_parameters["Number of points\ncalculated - nb_points_EIS"]["value"].get() < 0 ||
            typeof(choice_current_density_parameters["Power of the\ninitial frequency\n- f_power_min_EIS"]["value"].get()) != Int ||
            typeof(choice_current_density_parameters["Power of the\nfinal frequency\n- f_power_max_EIS"]["value"].get()) != Int ||
            typeof(choice_current_density_parameters["Number of frequencies\ntested - nb_f_EIS"]["value"].get()) != Int ||
            typeof(choice_current_density_parameters["Number of points\ncalculated - nb_points_EIS"]["value"].get()) != Int
        messagebox.showerror(title="f EIS", message="f_EIS parameters should be integer and number of points should be positive.")
        return
    end

    if choice_computing_parameters["Purge time - t_purge (s)"]["value"].get() < 0 ||
            choice_computing_parameters["Time between two purges\n- Δt_purge (s)"]["value"].get() < 0
        messagebox.showerror(title="Purge times", message="Negative times does not characterise purges.")
        return
    end

    if choice_computing_parameters["Number of GC nodes - nb_gc"]["value"].get() < 1 ||
            typeof(choice_computing_parameters["Number of GC nodes - nb_gc"]["value"].get()) != Int
        messagebox.showerror(title="nb_gc", message="The nb_gc value should be an integer bigger or equal to 1.")
        return
    end

    if choice_computing_parameters["Number of GDL nodes - nb_gdl"]["value"].get() < 1 ||
            typeof(choice_computing_parameters["Number of GDL nodes - nb_gdl"]["value"].get()) != Int
        messagebox.showerror(title="nb_gdl", message="The nb_gdl value should be an integer bigger or equal to 1.")
        return
    end

    if choice_computing_parameters["Number of MPL nodes - nb_mpl"]["value"].get() < 1 ||
            typeof(choice_computing_parameters["Number of MPL nodes - nb_mpl"]["value"].get()) != Int
        messagebox.showerror(title="nb_mpl", message="The nb_mpl value should be an integer bigger or equal to 1.")
        return
    end

    if choice_computing_parameters["Solver relative tolerance - rtol"]["value"].get() > 1e-3 ||
            choice_computing_parameters["Solver absolute tolerance - atol"]["value"].get() > 1e-3
        messagebox.showerror(title="Solver tolerance", message="rtol and atol should be lower than 1e-3 to limit the numerical errors.")
        return
    end

    if current_button == 0 && choice_buttons["type_display"]["value"].get() == 2 &&
            choice_buttons["type_plot"]["value"].get() == 1
        messagebox.showerror(title="n gdl", message="dynamic plot is not thought to be used with step current and multiple display. There would be too much plots to handle.")
        return
    end
end

"""
    set_equal_width(frame1::Any, frame2::Any, frame3::Any, frame4::Any, frame5::Any, frame6::Any)

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
function set_equal_width(frame1, frame2, frame3, frame4, frame5, frame6)
    # Initialisation of the list of widths
    widths = []

    for frame in [frame1, frame2, frame3, frame4, frame5, frame6]
        # Update the frame sizes
        frame.update_idletasks()
        # Get the current width of all frames
        push!(widths, frame.winfo_width())
    end

    # Set all frames to the maximum width
    for frame in [frame1, frame2, frame3, frame4, frame5, frame6]
        for i in 0:5
            frame.grid_columnconfigure(i, minsize=maximum(widths) / 5.5)  # Set minimum width of all column to max_width / 5.5
        end
    end
end

# ... Placeholder for launch_AlphaPEM_for_step_current, launch_AlphaPEM_for_polarization_current, launch_AlphaPEM_for_EIS_current
# These functions are large and would follow the same translation pattern as above.

# The following functions are commented out for now as they depend on the full AlphaPEM model implementation
# They will be fully translated when the model is available.

"""
    launch_AlphaPEM_for_step_current(simulator::Any, operating_inputs::Dict, current_parameters::Dict, computing_parameters::Dict)

Launch the AlphaPEM simulator for a step current density and display the results.

(Note: Full implementation would be added here - function signature is provided)
"""
# function launch_AlphaPEM_for_step_current(simulator, operating_inputs, current_parameters, computing_parameters)
#     # Starting time
#     start_time = _time()
#
#     # Figures preparation
#     fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)
#
#     # ... rest of implementation
# end

"""
    launch_AlphaPEM_for_polarization_current(simulator::Any, operating_inputs::Dict, current_parameters::Dict, computing_parameters::Dict)

Launch the AlphaPEM simulator for a polarization current density and display the results.

(Note: Full implementation would be added here - function signature is provided)
"""
# function launch_AlphaPEM_for_polarization_current(simulator, operating_inputs, current_parameters, computing_parameters)
#     # Starting time
#     start_time = _time()
#
#     # Figures preparation
#     fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)
#
#     # ... rest of implementation
# end

"""
    launch_AlphaPEM_for_EIS_current(simulator::Any, operating_inputs::Dict, current_parameters::Dict, computing_parameters::Dict)

Launch the AlphaPEM simulator for an EIS current density and display the results.

(Note: Full implementation would be added here - function signature is provided)
"""
# function launch_AlphaPEM_for_EIS_current(simulator, operating_inputs, current_parameters, computing_parameters)
#     # Starting time
#     start_time = _time()
#
#     # Figures preparation
#     fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)
#
#     # ... rest of implementation
# end

end  # end module GUIModules

