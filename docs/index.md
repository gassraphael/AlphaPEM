# AlphaPEM

AlphaPEM is an open-source software package for simulating proton exchange membrane fuel cell (PEMFC) systems for embedded applications. It is based on a physics-based, one-dimensional (1D), dynamic, two-phase, and non-isothermal model. It quickly simulates the internal states and voltage dynamics of PEMFC systems for all current densities and operating conditions imposed on it. In particular, it is possible to apply a step current density or use current profiles to generate polarization curves or electrochemical impedance spectroscopy (EIS) curves. It can also automatically calibrate the undetermined parameters of the model to simulate a given real fuel cell system.

A detailed [presentation](https://doi.org/10.48550/arXiv.2407.12373) of this program was published in the peer-reviewed journal SoftwareX. 

Improvements to **AlphaPEM**, such as the addition of heat transfer modelling and spatial extension to 1D+1D, will be available in the future.

Important note: AlphaPEM is an ongoing research project and is not a commercial product. Therefore, the latest online version may contain bugs, and not all features may be available. Relatively stable versions are listed in the [Major updates](#major-updates) section.


# Table of Contents

- [Installation](#installation)
- [Major updates](#major-updates)
- [Work in progress](#work-in-progress)
- [Related publications](#related-publications) 


# Installation

To install **AlphaPEM**, follow these steps in a shell:

1. Clone the repository:
    ```git clone https://github.com/gassraphael/AlphaPEM.git```

2. Navigate to the project directory:
    ```cd AlphaPEM```

3. Update the Python package manager, pip, to the latest available version:
    ```pip install --upgrade pip```

4. Install the required dependencies (eventually in a specific environment):
    ```
    pip install numpy scipy matplotlib colorama pygad
    python3 -m pip install git+https://github.com/RedFantom/ttkthemes
    ```
    

# Major updates

- V1.2 - in progress - This version of AlphaPEM includes: 
	- the addition of the MPL to the simulated cell, in both the anode and cathode. To ensure good numerical stability, a transition layer between each GDL and MPL is also included.
	- the addition of the open-source [ZSW GenStack](https://zenodo.org/records/14223364) as a calibrated fuel cell case study. 
- [V1.1](https://github.com/gassraphael/AlphaPEM/tree/11f07bd084a09cc6432f441b010d89d2a4229e4e) - 2025.08.18 - This version of AlphaPEM includes: 
	- the addition of heat transfer to the program, in cooperation with Pedro Affonso Nobrega (PERSEE, MINES PSL).
	- an improvement of the initial variable values: the algorithm waits for a given time (approximately 2 virtual hours) to reach equilibrium, and then the experiment starts (step/pola/EIS).
	- the limit liquid water saturation coefficient ($s_{lim}$) is temporarily removed for future refinement.
- [V1.0](https://github.com/gassraphael/AlphaPEM/tree/2b042c3d16d53fcd16779a5ffdc81eea75a9d012) - 2024.09.05 - This version of AlphaPEM corresponds to the one developed during Raphaël Gass's PhD from 2021 to 2024. 
	- It is based on a physics-based, one-dimensional (1D), dynamic, two-phase, and isothermal model.


# Work in progress

- The polarization curves from the EH-31 fuel cell example are no longer calibrated due to recent modifications made to the equations. A calibration of the undeterminate parameters will be performed in the future to correct this issue. If accurate examples are required, the [V1.0 version](https://github.com/gassraphael/AlphaPEM/tree/2b042c3d16d53fcd16779a5ffdc81eea75a9d012) of AlphaPEM can be used.


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
