# AlphaPEM Web Simulator

This directory contains the web-based graphical user interface for AlphaPEM, built entirely in Julia using Genie.jl and Stipple.jl.

## Structure

```
web/
├── app/                          # Genie application configuration
│   ├── app.jl                   # Server configuration
│   └── routes.jl                # HTTP routes and API endpoints
├── pages/                        # HTML templates
│   ├── index.html               # Landing page
│   ├── simulator.html           # Main simulator interface (tabs)
│   └── results.html             # Results display page
├── public/                       # Static assets
│   ├── css/
│   │   └── style.css            # Main stylesheet
│   └── js/                       # Client-side scripts (optional)
├── SimulatorBackend.jl          # Backend logic and AlphaPEM integration
├── WebApp.jl                    # Application entry point
└── README.md                    # This file
```

## Running the Simulator

### From command line (recommended for users):

```bash
cd /path/to/AlphaPEM
julia run_web_app.jl
```

The simulator will:
1. Start a local Genie server on `localhost:8000`
2. Automatically open your default browser
3. Display the AlphaPEM interface

### From Julia REPL:

```julia
using AlphaPEM
AlphaPEM.Interfaces.launch_web_simulator()
```

## File Descriptions

### `app/app.jl`
Configures the Genie web server:
- Server host (localhost only for security)
- Port number (8000)
- Worker threads
- Logging level

### `app/routes.jl`
Defines HTTP endpoints:
- `GET /` - Home page
- `GET /simulator` - Simulator interface
- `GET /results/:result_id` - Results page
- `POST /api/simulate/step` - Launch step simulation
- `POST /api/simulate/polarization` - Launch polarization simulation
- `POST /api/simulate/eis` - Launch EIS simulation
- `GET /api/status/:result_id` - Check simulation status
- `GET /api/results/:result_id` - Retrieve results data

### `SimulatorBackend.jl`
Backend module providing:
- `get_available_fuel_cells()` - List available models
- `get_fuel_cell_defaults(fuel_cell_type)` - Default parameters
- `validate_parameters(params)` - Input validation
- `build_simulation_config(params, sim_type)` - Create SimulationConfig
- `run_step_simulation(config)` - Execute step simulation
- `run_polarization_simulation(config)` - Execute polarization curve
- `run_eis_simulation(config)` - Execute EIS simulation
- `get_simulation_status(result_id)` - Check status
- `get_simulation_results(result_id)` - Retrieve results

### `WebApp.jl`
Main entry point:
- Loads configuration
- Initializes backend
- Loads routes
- Starts Genie server

### HTML Pages

#### `index.html`
Landing page with:
- Project overview
- Available simulation types
- Links to documentation
- Call-to-action button to simulator

#### `simulator.html` (Development in Étape 2-3)
Main interface with:
- Tabbed navigation:
  - Fuel Cell Selection
  - Operating Conditions
  - Accessible Parameters
  - Undetermined Parameters
  - Model Configuration
  - Computing Parameters
  - Simulation Selection & Launch
- Form validation
- Parameter presets

#### `results.html` (Development in Étape 4)
Results display with:
- Simulation summary
- Plots visualization
- Data export options
- Links to new simulations

### `public/css/style.css`
Comprehensive stylesheet providing:
- Responsive design (mobile-friendly)
- Color scheme and typography
- Component styles (buttons, forms, cards)
- Utilities and helper classes
- Dark mode support (optional)

## API Endpoints (Detailed)

### GET `/api/fuel_cells`
Returns available fuel cell models:
```json
{
  "default": {"name": "Default Model", "description": "..."},
  "ZSW_GenStack": {...},
  ...
}
```

### GET `/api/fuel_cell_defaults/:fuel_cell_type`
Returns default parameters for fuel cell type.

### POST `/api/validate_parameters`
Validates parameters before simulation.

### POST `/api/simulate/step`
Launches step current simulation.

### POST `/api/simulate/polarization`
Launches polarization curve simulation.

### POST `/api/simulate/eis`
Launches EIS spectroscopy simulation.

### GET `/api/status/:result_id`
Returns current simulation status (running/completed/failed).

### GET `/api/results/:result_id`
Returns simulation results data.

## Security Considerations

- Server binds to `127.0.0.1:8000` (localhost only)
- No automatic network exposure
- Parameters validated on both client and server
- Results stored locally in `/results/`
- No data sent to external servers

## Logging

Genie logs are displayed in console:
- Log level: "warn" in production, "info" in development
- Logs include route hits, errors, simulation progress

