# AlphaPEM

AlphaPEM is an open-source software package for simulating proton exchange membrane fuel cell (PEMFC) systems for embedded applications. It is based on a physics-based, one-dimensional (1D), dynamic, two-phase, and isothermal model. It can quickly simulate the internal states and voltage dynamics of PEMFC systems, and produce polarization and EIS curves. It can also automatically calibrate the undetermined parameters of the model to simulate a given real fuel cell system.

A detailed [presentation](https://doi.org/10.48550/arXiv.2407.12373) of this program was published in the peer-reviewed journal SoftwareX. 

Improvements to **AlphaPEM**, such as the addition of heat transfer modelling and spatial extension to 1D+1D, will be available in the future.


# Table of Contents

- [Installation](#installation)
- [Major updates](#major-updates)
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
    pip install numpy scipy matplotlib colorama geneticalgorithm2
    python3 -m pip install git+https://github.com/RedFantom/ttkthemes
    ```
    

# Major updates

- V1.0 - 2024.09.05 - This version of AlphaPEM corresponds to the one developed during Raphaël Gass's PhD from 2021 to 2024.


# Related publications

The detailed model description and simulation results can be found in the following articles and thesis.

- Submitted journal papers:
    - **AlphaPEM: An Open-Source Dynamic 1D Physics-Based Pem Fuel Cell Model for Embedded Applications**
        - In the journal SoftwareX, in [arXiv](https://doi.org/10.48550/arXiv.2407.12373), in [HAL](https://hal.science/hal-04647829) or in [SSRN](http://ssrn.com/abstract=4946674) (preprint).
        - The objective of this work is to highlight the AlphaPEM software, which has been published as open-source on GitHub. The first version of this PEM fuel cell simulator is based on the dynamic 1D model developed during 2021-2024.

	- **An Advanced 1D Physics-Based Model for PEM Hydrogen Fuel Cells With Enhanced Overvoltage Prediction**
		- In the International Journal of Hydrogen Energy, in [arXiv](https://doi.org/10.48550/arXiv.2404.07508) or in [HAL](https://hal.science/hal-04530852) (preprint).
		- The aim of this study was to introduce the dynamic 1D model developed during 2021-2024, emphasizing the adjustment of the equations for this specific model and their numerical resolution. Furthermore, a novel coefficient is proposed to establish a physical relationship between the voltage drop at high currents, the quantity of liquid water in the cell, and operating conditions.
	
- Published journal papers:
	- **A Critical Review of Proton Exchange Membrane Fuel Cells Matter Transports and Voltage Polarisation for Modelling**
		- In the [Journal of the Electrochemical Society](https://doi.org/10.1149/1945-7111/ad305a) or in [HAL](https://hal.science/hal-04493419) (postprint).
		- The aim of this work was to compile, within a single article, all the equations required for the physical modeling of a fuel cell. Each equation is complemented with explanations, critical analysis, and suggestions for potential enhancements.
		
- Thesis:
	- to complete.
