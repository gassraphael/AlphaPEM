# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from modules.main_modules import launch_AlphaPEM_for_step_current, launch_AlphaPEM_for_polarization_current, \
                                 launch_AlphaPEM_for_polarization_current_for_calibration, \
                                 launch_AlphaPEM_for_EIS_current

# __________________________________________________AlphaPEM settings___________________________________________________
"""
- Users can select various preconfigured configurations for execution.
- Adjustments to these configurations can be made within setting.py or current_densities.py and their associated files.
- Selecting different type_fuel_cell during a single run results in simultaneous plots for various configurations.
"""

if __name__ == '__main__':
    # Fuel cell possibilities: "EH-31_1.5"(2021), "EH-31_2.0"(2021), "EH-31_2.25"(2021), "EH-31_2.5"(2021),
    #                          "LF"(2010), "manual_setup".
    # This parameter includes the fuel cell used in the model and the corresponding operating conditions.
    # - EH-31 is a fuel cell developed by EH GROUP. 1.5, 2.0, 2.25 and 2.5 corresponds to the different pressure options.
    # - LF corresponds to the fuel cell used in Linhao Fan work: http://dx.doi.org/10.1016/j.enconman.2017.08.034.
    type_fuel_cell_1 = "EH-31_2.0"
    type_fuel_cell_2 = None
    type_fuel_cell_3 = None
    type_fuel_cell_4 = None
    # Current density possibilities: "step", "polarization", "polarization_for_cali", "EIS".
    type_current = "polarization"
    # Auxiliary system possibilities: "forced-convective_cathode_with_anodic_recirculation",
    #                                 "forced-convective_cathode_with_flow-through_anode", "no_auxiliary".
    type_auxiliary = "forced-convective_cathode_with_flow-through_anode"
    # Control strategy for the operating inputs: "Phi_des", "no_control".
    type_control_1 = "no_control"
    type_control_2 = "no_control"
    type_control_3 = "no_control"
    type_control_4 = "no_control"
    # Purges possibilities: "constant_purge", "periodic_purge", "no_purge".
    type_purge = "no_purge"
    # Display possibilities: "multiple", "synthetic", "no_display".
    type_display = "multiple"
    # Plot possibilities: "dynamic", "fixed". Using dynamic plot option enables real-time figure updates during program
    # execution, albeit at the cost of decreased program speed.
    type_plot = "fixed"

# __________________________________________________________Main________________________________________________________
    """This section is dedicated to ensuring the proper execution of the simulator, considering all the various 
    possibilities including real-time figure updates and simultaneous plotting for different configurations. 
    This should remain unaltered for regular program usage.
    """

    # Check if the type_current is valid
    if type_current == "step":
        launch_AlphaPEM_for_step_current(type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3, type_fuel_cell_4,
                                         type_current, type_auxiliary, type_control_1, type_control_2, type_control_3,
                                         type_control_4, type_purge, type_display, type_plot)
    elif type_current == "polarization":
        launch_AlphaPEM_for_polarization_current(type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3, type_fuel_cell_4,
                                                 type_current, type_auxiliary, type_control_1, type_control_2,
                                                 type_control_3, type_control_4, type_purge, type_display, type_plot)
    elif type_current == "polarization_for_cali":
        launch_AlphaPEM_for_polarization_current_for_calibration(type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3,
                                                                 type_fuel_cell_4, type_current, type_auxiliary,
                                                                 type_control_1, type_control_2, type_control_3,
                                                                 type_control_4, type_purge, type_display, type_plot)
    elif type_current == "EIS":
        launch_AlphaPEM_for_EIS_current(type_fuel_cell_1, type_current, type_auxiliary, type_control_1, type_purge,
                                        type_display, type_plot)
    else:
        raise ValueError('You have to specify a type_current which is accepted.')

