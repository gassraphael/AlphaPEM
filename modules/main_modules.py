# -*- coding: utf-8 -*-

"""This module contains some of the required functions for the main.py file.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import os
import matplotlib as mpl
import matplotlib.pyplot as plt


# _____________________________________________________Main modules_____________________________________________________

def figures_preparation(type_current, type_display):
    """ This function create the required figures and axes according to the type_current and type_display.

    Parameters
    ----------
    type_current : str
        Type of current density function.
    type_display : str
        Type of display.

    Returns
    -------
    fig1 : matplotlib.figure.Figure
        Figure for the first plot.
    ax1 : matplotlib.axes._subplots.AxesSubplot
        Axes for the first plot.
    fig2 : matplotlib.figure.Figure
        Figure for the second plot.
    ax2 : matplotlib.axes._subplots.AxesSubplot
        Axes for the second plot.
    """

    mpl.rcParams['font.family'] = 'cmr10'  # 'cmr10' for English characters and 'DejaVu Serif' for French ones
    mpl.rcParams['axes.formatter.use_mathtext'] = True  # For the scientific notation
    mpl.rcParams['lines.linewidth'] = 2.0
    mpl.rcParams['lines.markersize'] = 5.0

    if type_display == "no_display":
        fig1, ax1 = None, None
        fig2, ax2 = None, None
        fig3, ax3 = None, None

    # For the step current
    if type_current == "step":
        if type_display == "multiple":  # saving instruction is directly implemented within AlphaPEM.Display here.
            mpl.rcParams['font.size'] = 18  # Font size for all text
            fig1, ax1 = None, None  # Here, additional plots are unnecessary
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
        elif type_display == "synthetic":
            mpl.rcParams['font.size'] = 13  # Font size for all text
            fig1, ax1 = plt.subplots(3, 3, figsize=(14, 14))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
            plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)

    # For the polarization curve
    elif type_current == "polarization":
        if type_display == "multiple":
            mpl.rcParams['font.size'] = 11  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3, figsize=(14, 4.7))
            fig2, ax2 = plt.subplots(1, 3, figsize=(14, 4.7))
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
            plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        elif type_display == "synthetic":
            mpl.rcParams['font.size'] = 18  # Font size for all text
            fig1, ax1 = plt.subplots(figsize=(8, 8))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary

    # For the EIS curve
    elif type_current == "EIS":
        if type_display == "multiple":
            mpl.rcParams['font.size'] = 18  # Font size for all text
            fig1, ax1 = plt.subplots(figsize=(8, 8))
            fig2, ax2 = plt.subplots(figsize=(8, 8))
            fig3, ax3 = plt.subplots(figsize=(8, 8))
        elif type_display == "synthetic":
            mpl.rcParams['font.size'] = 13  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3, figsize=(14, 4.7))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
            plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)

    return fig1, ax1, fig2, ax2, fig3, ax3


def plot_saving(type_fuel_cell, type_current, type_display, fig1, fig2, fig3):
    """This function saves the plots in the mentioned folder. The names of the files are automatically generated
    according to the type_current and the type_display.

    Parameters
    ----------
    type_fuel_cell : str
        Type of the fuel cell.
    type_current : str
        Type of current density function.
    type_display : str
        Type of display.
    fig1 : matplotlib.figure.Figure
        Figure for the first plot.
    fig2 : matplotlib.figure.Figure
        Figure for the second plot.
    """

    # Folder name
    subfolder_name = type_fuel_cell[:type_fuel_cell.rfind('_')] if type_fuel_cell.rfind('_') != -1 else type_fuel_cell

    # For the step current
    if type_current == "step":
        if type_display == "multiple":
            pass  # saving instruction is directly implemented within AlphaPEM.Display for this situation.
        if type_display == "synthetic":
            saving_instructions("results", subfolder_name, "step_current_syn_1.pdf", fig1)

    # For the polarization curve
    elif type_current == "polarization":
        if type_display == "multiple":
            saving_instructions("results", subfolder_name, "global_indicators_1.pdf", fig1)
            saving_instructions("results", subfolder_name, "pola_curve_syn_1.pdf", fig2)
        elif type_display == "synthetic":
            saving_instructions("results", subfolder_name, "pola_curve_1.pdf", fig1)

    # For the EIS curve
    elif type_current == "EIS":
        if type_display == "multiple":
            saving_instructions("results", subfolder_name, "Nyquist_plot_1.pdf", fig1)
            saving_instructions("results", subfolder_name, "Bode_amplitude_curve_1.pdf", fig2)
            saving_instructions("results", subfolder_name, "Bode_angle_curve_1.pdf", fig3)
        elif type_display == "synthetic":
            saving_instructions("results", subfolder_name, "Nyquist_plot_syn_1.pdf", fig1)


def saving_instructions(root_folder, subfolder_name, filename, fig):
    """This function gives the saving instructions for the figures.

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
