# AlphaPEM

*Advanced Fuel Cell Model for Embedded Applications*

## Overview

AlphaPEM is an open-source Julia package for simulating proton exchange membrane fuel cell (PEMFC) systems in 
embedded and control applications. Built on a physics-based, finite-volume, pseudo-two-dimensional (1D+1D) dynamic 
model, it accurately captures internal states and voltage dynamics across all operating conditions.

**Key Capabilities:**
- **Physics-Based Modeling**: Finite-volume pseudo-2D model with two-phase and non-isothermal dynamics.
- **Flexible Configuration**: Support for co-flow and counter-flow gas channels, multiple fuel cell models.
- **Characterization**: Generate polarization curves and electrochemical impedance spectroscopy (EIS) measurements.
- **Automatic Calibration**: Genetic algorithm-based parameter identification for real fuel cell systems.
- **Modern Interface**: Web-based graphical interface (V2.0+) or programmatic access via command line interface (CLI).

## Quick Start

### Installation (from source)

```bash
git clone https://github.com/gassraphael/AlphaPEM.git && cd AlphaPEM
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run a simulation

**Via Web Interface:**
```bash
julia --project=. run_web_app.jl
```
Then navigate to `http://localhost:8000` in your browser.

**Via Command Line:**
```bash
julia --project=. examples/run_step.jl
```

## Documentation sections

- **[Getting Started](getting_started/installation.md)** — Complete installation guide and initial setup.
- **[User Guide](user_guide/web_interface.md)** — How to use the web interface, CLI examples, and workflows.
- **[Advanced Topics](advanced/calibration.md)** — Parameter calibration and validity analysis.
- **[About](about/roadmap.md)** — Project roadmap and publications.

## Features by version

**Version 2.0** (in development)
- Transition to Julia for improved performance, maintainability, and modern object-oriented architecture.
- New web-based interface replacing the previous GUI.
- Enhanced O₂-to-Pt particle flow for high-current overvoltage modeling.
- Spatial extension to 1D+1D with convective flow in GCs (thermal evolution remains 1D), including liquid water transport in gas channels with sorption flow.
- Added heat transfer (with Pedro Affonso Nobrega, PERSEE, MINES PSL) and MPL integration in both anode and cathode.

**[Version 1.0](https://github.com/gassraphael/AlphaPEM/tree/2b042c3d16d53fcd16779a5ffdc81eea75a9d012)** (2024-09-05)
- Initial release: 1D dynamic, two-phase, isothermal PEMFC model.

## Peer-reviewed publication

The AlphaPEM V1.0 model and software are described in detail in:

> **AlphaPEM: An Open-Source Dynamic 1D Physics-Based PEM Fuel Cell Model for Embedded Applications**  
> Raphaël Gass, 2025  
> *SoftwareX* | [DOI](https://doi.org/10.1016/j.softx.2024.102002) | [HAL](https://hal.science/hal-04647829)

Additional physics and methodology papers are available in the [Publications](about/publications.md) section.

## For Developers

- Use AlphaPEM as a Julia package in your projects with full API access.
- Integrate PEMFC simulations into your Python applications via `juliacall`.
- View the [CLI Usage Guide](user_guide/cli_usage.md) for programmatic access and integration examples.

## Support and Contact

For questions, issues, or collaboration inquiries, please contact:
📧 **raphael.gass@univ-reunion.fr**

**Important Note**: AlphaPEM is an active research project. Newer versions may be available pending scientific 
publication; contact us to discuss collaboration opportunities.

---

**License**: [GNU GPL 3.0](https://github.com/gassraphael/AlphaPEM/blob/main/LICENSE)
