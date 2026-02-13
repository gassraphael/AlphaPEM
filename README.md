# AlphaPEM

AlphaPEM is an open-source software package for simulating proton exchange membrane fuel cell (PEMFC) systems for embedded applications. It is based on a physics-based, finite-volume, pseudo-two-dimensional (1D+1D), dynamic, two-phase, and non-isothermal model. It quickly simulates the internal states and voltage dynamics of PEMFC systems for all current densities and operating conditions imposed on it. In particular, it is possible to apply a step current density or use current profiles to generate polarization curves or electrochemical impedance spectroscopy (EIS) curves. It can also automatically calibrate the undetermined parameters of the model to simulate a given real fuel cell system.

A detailed [presentation](https://doi.org/10.48550/arXiv.2407.12373) of this program has been published in the peer-reviewed journal SoftwareX (limited to [version V1.0](#major-updates)). Furthermore, comprehensive [documentation](https://gassraphael.github.io/AlphaPEM/) in Numpy style for the software functions is available.

Improvements to AlphaPEM are discussed in the [roadmap section](#roadmap).

**Important note:** AlphaPEM is an ongoing research project and is not a commercial product. Therefore, the latest online version may contain bugs, and not all features may be available. The current work is detailled in the [work in progress](#work-in-progress) section. Relatively stable versions are listed in the [Major updates](#major-updates) section.

![AlphaPEM graphical user interface](docs/images/demo.png "AlphaPEM graphical user interface (GUI)")

# Table of Contents

- [Installation](#installation)
- [Start](#start)
- [Major updates](#major-updates)
- [Work in progress](#work-in-progress)
- [Roadmap](#roadmap)
- [Related publications](#related-publications) 
- [Contributions](#contributions)
- [Contact](#contact)


# Installation

To install **AlphaPEM**, follow these steps in a shell:

1. Clone the repository:
    ```sh
    git clone https://github.com/gassraphael/AlphaPEM.git
    ```

2. Navigate to the project directory:
    ```sh
    cd AlphaPEM
    ```
    
3. Update the Python package manager, pip, to the latest available version:
    ```sh
    pip install --upgrade pip
    ```
    
4. Create a new environment, referred to here as *env*, and activate it:
    ```sh
    python3 -m venv env
    source env/bin/activate
    ```

4. Install the required dependencies (eventually in a specific environment):
    ```sh
    pip install numpy scipy matplotlib colorama pygad
    python3 -m pip install git+https://github.com/RedFantom/ttkthemes
    ```

# Start

You have two main ways to run AlphaPEM:

## Using the graphical user interface (GUI)

The GUI provides a quick way to configure and run simulations without modifying the source code. However, it does not yet grant access to all the functionalities of the code.

1. Execute the GUI file:

   ```sh
   python3 GUI.py
   ```

2. In the GUI (as shown in the figure of [AlphaPEM section](#alphapem)):

   - Select a predefined fuel cell specification from the 'Fuel cell' dropdown menu. Operating conditions and parameters can also be adjusted by selecting 'Enter your specifications' in this menu.

   - Choose the configuration you prefer for the auxiliaries, voltage zone, purge, display and plot under 'Model configuration'.

   - Select your desired simulation type at the bottom of the GUI (e.g., current density step, polarization curve, or EIS curve).

3. Run the simulation to generate results (internal states and voltage dynamics) in the /results directory.

## Using the command line (programmers)

The main.py file is used for standard operation and provides full control for programmers. This allows for using any physically acceptable current density function, beyond the predefined configurations of the GUI.

1. Modify parameters and input current densities directly in the appropriate configuration files (e.g., /configuration/settings.py or /configuration/current_densities.py) if needed.

2. Select a predefined fuel cell specification, a given configuration and the desired simulation directly at the beguining of the /main.py file.

3. Execute the main file to generate results (internal states and voltage dynamics) in the /results directory:

```sh
python3 main.py
```

4. Automated parameter calibration (advanced)

   To adapt AlphaPEM to a new, specific fuel cell, you must calibrate the undetermined physical parameters (like GDL porosity) using experimental data. This functionality is not yet available from the GUI. The calibration uses a genetic algorithm (PyGAD) to match simulated results to experimental data.

   1. Input Experimental Data: place experimental polarization curves (at least three) into the file: ./calibration/experimental_values.

   2. Configure Parameters: input the operating conditions and accessible physical parameters of your fuel cell system in: ./modules/calibration_modules.

   3. Run Calibration: execute the calibration program (preferably on a computing cluster due to computational cost):
   ```sh
   python3 ./calibration/parameter_calibration.py
   ```

# Major updates

- V1.3 - in progress - This version of AlphaPEM includes: 
	- the addition of O2 flow to Pt particules which improves the modeling of overvoltage due to flooding at high curent densities.
		- the limiting liquid water saturation coefficient ($s_{lim}$) has been definitively removed, as this model replaces it.
	- the addition of liquid water flow inside the GC (with the sorption flow at the GDL/GC interface).
	- the spatial extension to 1D+1D (except thermal evolution which remains 1D for now).
- [V1.2](https://github.com/gassraphael/AlphaPEM/tree/b71f42878a186e17efeb7e97b5d7fb50d6e76827) - 2025.12.11 - This version of AlphaPEM includes: 
	- the addition of convective flow between the inlet, gas channels, and outlet of the cell, thereby removing the Pukrushpan equations (from Michigan University).
		- auxiliaries are temporarily removed, as they require reconstruction. 
	- the addition of the MPL to the simulated cell, in both the anode and cathode. 
	- effective diffusive flows for the dissolved water insided the CLs are introduced.
	- the addition of the open-source [ZSW GenStack](https://zenodo.org/records/14223364) as a calibrated fuel cell case study. 
- [V1.1](https://github.com/gassraphael/AlphaPEM/tree/11f07bd084a09cc6432f441b010d89d2a4229e4e) - 2025.08.18 - This version of AlphaPEM includes: 
	- the addition of heat transfer to the program, in cooperation with Pedro Affonso Nobrega (PERSEE, MINES PSL).
	- an improvement of the initial variable values: the algorithm waits for a given time to reach equilibrium, and then the experiment starts (step/pola/EIS).
	- the limiting liquid water saturation coefficient ($s_{lim}$) is temporarily removed for future refinement.
- [V1.0](https://github.com/gassraphael/AlphaPEM/tree/2b042c3d16d53fcd16779a5ffdc81eea75a9d012) - 2024.09.05 - This version of AlphaPEM corresponds to the one developed during Raphaël Gass's PhD from 2021 to 2024. 
	- It is based on a physics-based, one-dimensional (1D), dynamic, two-phase, and isothermal model.


# Work in progress

- Calibration of the model using pre-selected data from ZSW-GenStack or EH-31 is currently underway.
- Auxiliaries are temporarily removed, as they require reconstruction. 


# Roadmap

- Spatial extension to 1D+1D for modeling the thermal evolution.
- Spatial extension to 1D+1D+1D: a 1D channel will be added to each manifold, enabling full-stack modeling.
- Integration of more accurate physical models for the auxiliaries.
- Inclusion of ECSA degradation in the simulation framework.
- Enhancement of the GUI to allow seamless addition of new fuel cell configurations without modifying the source code.

# Related publications

The detailed model description and simulation results can be found in the following articles and thesis.
	
- Published journal papers:
	- **AlphaPEM: An Open-Source Dynamic 1D Physics-Based Pem Fuel Cell Model for Embedded Applications** (2025, 1st author)
	    - In the [SoftwareX](https://doi.org/10.1016/j.softx.2024.102002) journal, in [arXiv](https://doi.org/10.48550/arXiv.2407.12373), in [HAL](https://hal.science/hal-04647829) or in [SSRN](http://ssrn.com/abstract=4946674) (postprint).
	    - The objective of this work is to highlight the AlphaPEM software, which has been published as open-source on GitHub. The first version of this PEM fuel cell simulator is based on the dynamic 1D model developed during 2021-2024. 

	- **An Advanced 1D Physics-Based Model for PEM Hydrogen Fuel Cells With Enhanced Overvoltage Prediction** (2025, 1st author)
		- In the [International Journal of Hydrogen Energy](https://doi.org/10.1016/j.ijhydene.2024.11.374), in [arXiv](https://doi.org/10.48550/arXiv.2404.07508), in [HAL](https://hal.science/hal-04530852) or in [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4812343) (postprint).
		- The aim of this study was to introduce the dynamic 1D model developed during 2021-2024, emphasizing the adjustment of the equations for this specific model and their numerical resolution. Furthermore, a novel coefficient is proposed to establish a physical relationship between the voltage drop at high currents, the quantity of liquid water in the cell, and operating conditions.
		- ![1D modeling of matter transport phenomena in a PEM single cell divided into several nodes.](docs/images/nodal_model.png "1D modeling of matter transport phenomena in a PEM single cell divided into several nodes")
		
	
	- **A Critical Review of Proton Exchange Membrane Fuel Cells Matter Transports and Voltage Polarisation for Modelling** (2024, 1st author)
		- In the [Journal of the Electrochemical Society](https://doi.org/10.1149/1945-7111/ad305a) or in [HAL](https://hal.science/hal-04493419) (postprint).
		- The aim of this work was to compile, within a single article, all the equations required for the physical modeling of a fuel cell. Each equation is complemented with explanations, critical analysis, and suggestions for potential enhancements.
		
- Thesis manuscript:
	- **Advanced physical modeling of PEM fuel cells to enhance their performances** (2024, 1st author)
		- In [HAL](https://hal.science/tel-04923016) (final version).
		- The objective of this thesis was to develop an advanced model for PEMFCs to optimize their control and improve performance. A 1D, dynamic, two-phase, isothermal model was proposed, leading to the development of the open-source software AlphaPEM, which enables accurate simulations and facilitates predictive control strategies for enhanced fuel cell operation.


# Contributions

## Authors

- AlphaPEM is firstly developed by [Raphaël Gass](https://gassraphael.github.io/) during his PhD thesis in control engineering at the [LIS Laboratory](https://www.lis-lab.fr/) in [Aix-Marseille University](https://www.univ-amu.fr/), and in co-supervision with [FEMTO-ST Institute](https://www.femto-st.fr/en), within the [FCLab](https://www.fclab.fr/), in [Franche-Comté University](https://www.univ-fcomte.fr/), from 2021 to 2024. This work has been supervised by Prof. Zhongliang Li (FEMTO-ST), Prof. Rachid Outbib (LIS), Prof. Samir Jemei (FEMTO-ST) and Prof. Daniel Hissel (FEMTO-ST).

- The development of AlphaPEM was subsequently continued by Raphaël Gass during his postdoctoral research from 2025 to 2027 at [ENERGY-Lab](https://www.energylab.re/), [University of Reunion island](https://www.univ-reunion.fr/), in partnership with the [ZSW Institute](https://www.zsw-bw.de/) in Ulm, Germany. This work was supervised by Prof. Michel Benne (ENERGY-Lab), Associate Prof. Cédric Damour (ENERGY-Lab), Associate Prof. Dominique Grondin (ENERGY-Lab), and Dr. Florian Wilhelm (ZSW).
    
## Financial support

This work has been supported:

- from 2021 to 2024 by French National Research Agency via project [DEAL](https://deal.lis-lab.fr/) (Grant no. ANR-20-CE05-0016-01), the Region Provence-Alpes-Côte d’Azur, the EIPHI Graduate School (contract ANR-17-EURE-0002) and the Region Bourgogne Franche-Comté.
- from 2025 to 2027 by European FEDER funds via project [OPUS-H2](https://www.energylab.re/projets/projets-en-cours/opus-h2/) and the Region Reunion.

## Licenses

**AlphaPEM** is licensed under the GNU GPL 3.0. See the [LICENSE](LICENSE) file for more details. 

It also includes components licensed under the [BSD-3-Clause license](calibration/LICENSE-BSD-3-CLAUSE):

- calibration/parameter_calibration.py from [PyGAD](https://github.com/ahmedfgad/GeneticAlgorithmPython). 

## New contributors

Contributions from the community are welcome! If you would like to contribute to **AlphaPEM**, please follow these steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature/YourFeature`).
3. Commit your changes (`git commit -am 'Add some feature'`).
4. Push to the branch (`git push origin feature/YourFeature`).
5. Create a new Pull Request.


# Contact

For any questions or support, please contact me at [gassraphael@proton.me](mailto:gassraphael@proton.me).

Thank you for using **AlphaPEM**!

