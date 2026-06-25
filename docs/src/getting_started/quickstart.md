# Quick start

Get AlphaPEM running in 5 minutes with these essential steps.

## 1. Install and run

```bash
# Clone the repository
git clone https://github.com/gassraphael/AlphaPEM.git && cd AlphaPEM

# Instantiate Julia environment (one-time setup)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Launch the web interface
julia --project=. run_web_app.jl
```

The web interface will automatically open at `http://localhost:8000`.

## 2. Run your first simulation (web interface)

1. **Home Tab**: Review AlphaPEM capabilities.
2. **Simulator Tab**:
   - **Model Selection**: Choose your stack (e.g., `ZSW_GenStack`, `EH31`, or `custom`).
   - **Physical parameters**: Set the system parameters, such as membrane thickness, electrode properties, and flow field design.
   - **Operating Conditions**: Set the temperature, pressures, stoichiometries, and humidities.
   - **Numerical parameters**: Adjust the number of nodes (allows to switch between 1D and 1D+1D simulations) and solver tolerances.
   - **Configuration**: Choose between different configurations for auxiliary equipment and purge strategies.
   - **Simulation Type**: 
     - **Step**: Temporal response to step current changes.
     - **Polarization**: Steady-state V-I curve.
     - **EIS**: Impedance spectroscopy (Nyquist/Bode).
3. **Results Tab**: Visualize dynamic plots and click **"Download CSV"** to export raw data.


<img src="docs/src/images/interface.png" alt="AlphaPEM web-based interface" width="100%" />

## 3. Run your first simulation (command line)

For more control, use example scripts:

```bash
# Step response at constant current density
julia --project=. examples/run_step.jl

# Polarization curve (voltage vs. current)
julia --project=. examples/run_polarization.jl

# Electrochemical impedance spectroscopy
julia --project=. examples/run_EIS.jl
```

Results are saved to the `results/` directory.

## 4. Customize a simulation

Edit an example script to change parameters:

```bash
# Open the step response example
nano examples/run_step.jl
```

Key configuration fields:

| Parameter | Options | Purpose |
|-----------|---------|---------|
| `type_fuel_cell` | `:ZSW_GenStack`, `:EH31`, etc. | Select fuel cell model |
| `type_current` | `StepParams(...)`, `PolarizationParams(...)`, `EISParams(...)` | Define current profile |
| `voltage_zone` | `:before_voltage_drop`, `:full` | Include voltage drop or not |
| `type_auxiliary` | `:no_auxiliary`, `:forced_convective_cathode_...` | Auxiliary equipment |
| `type_purge` | `:no_purge`, `:constant_purge`, `:periodic_purge` | Purge strategy |

## 5. Next steps

- **Full Installation**: See [Installation Guide](installation.md) for optional features
- **Web Interface**: Explore the [Web Interface Guide](../user_guide/web_interface.md)
- **Advanced Topics**: Learn about [Parameter Calibration](../advanced/calibration.md)

## Troubleshooting quick reference

| Problem | Solution |
|---------|----------|
| Slow first run | Julia precompilation is normal; subsequent runs are fast |
| Import error | Run `julia --project=. -e 'using Pkg; Pkg.update()'` |

---

**Ready to dive deeper?** Check out the [User Guide](../user_guide/web_interface.md) or contact us by [e-mail](mailto:raphael.gass@univ-reunion.fr).
