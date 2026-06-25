# Web interface guide

The AlphaPEM web interface provides an intuitive, browser-based environment for configuring and running PEMFC 
simulations without writing code. This guide covers all features and workflows.

## Launching the interface

Start the web server:

```bash
julia --project=. run_web_app.jl
```

The interface automatically opens in your default browser at `http://localhost:8000`. If not, navigate manually.

## Interface overview

The web interface is organized into four main sections:

### Home tab

Displays an overview of AlphaPEM's capabilities:

- Brief description of the PEMFC model
- Key features and simulation types
- Links to documentation and external resources
- General project information

### Simulator tab

Configure fuel cell model and operating conditions.

**Fuel cell model selection:**
- Dropdown menu with predefined models:
  - `:ZSW_GenStack` — Calibrated ZSW GenStack (open-source hardware stack)
  - `:EH31` — Calibrated EH31 fuel cell model
  - Custom configurations

**Physical and undetermined parameters:**

- Cell geometric parameters: such as active area, channel dimensions, and membrane thickness.
- Material properties: such as gas diffusion layer porosity, and interfacial resistance coefficient of O2 adsorption on the Pt sites.
- Electrochemical parameters: such as reference exchange current densities, and crossover correction coefficient.
- Manifold geometric parameters: such as inlet and exhaust throttle areas, and manifold length.


**Operating conditions:**

| Parameter | Unit | Description |
|-----------|------|-------------|
| Temperature | °C   | Cell temperature |
| Pressure (Anode) | bar  | Hydrogen inlet pressure |
| Pressure (Cathode) | bar  | Air inlet pressure |
| Stoichiometry (Anode) | —    | H₂ excess ratio |
| Stoichiometry (Cathode) | —    | Air excess ratio |
| Relative Humidity | %    | Inlet gas humidity |

**Numerical parameters:**

- **Mesh nodes**: Number of spatial nodes, which allows switching from 1D to 1D+1D simulations.
- **Solver tolerances**: Absolute and relative tolerances for DAE solver

**Model configuration:**

| Setting             | Options                                                      | Effect                             |
|---------------------|--------------------------------------------------------------|------------------------------------|
| Voltage zone        | `before voltage drop`, `full`                                | Include or exclude ohmic drops     |
| Auxiliary equipment | `no auxiliary`, `anodic recirculation`, `flow-through anode` | On-board or laboratory application |
| Flow type           | `counter flow`, `co-flow`                                    | Thermal and matter management      |
| Purge strategy      | `no_purge`, `constant_purge`, `periodic_purge`               | H2 management                      |

### Simulation tab

Configure and launch simulations.

**Simulation type selection:**

1. **Step response**
   - Apply a step current density
   - Observe system transient response

2. **Polarization curve**
   - Sweep current from low to high
   - Characterize fuel cell performance

3. **Electrochemical impedance spectroscopy (EIS)**
   - Apply small-signal AC perturbations
   - Identify internal impedance sources

### Results tab

View and export simulation outputs.

**Data export:**

- **PDF format**: export the current plot. 
- **XLSX format**: export all simulation data to Excel.

## Next steps

- **Batch simulations**: Use [Command Line Interface](cli_usage.md) for automated workflows.
- [Valid Parameter Region Analysis](validity_analysis.md) — Identify valid parameter space before calibration
- **Parameter calibration**: See [Calibration Guide](../advanced/calibration.md).

---

**Questions?** Contact [raphael.gass@univ-reunion.fr](mailto:raphael.gass@univ-reunion.fr) or open an issue on [GitHub](https://github.com/gassraphael/AlphaPEM/issues).
