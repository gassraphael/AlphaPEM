# Installation

This guide covers all installation methods for AlphaPEM, from basic development setup to advanced optional features.

## System Requirements

- **Julia**: 1.11 or later (download from [julialang.org](https://julialang.org/downloads/))
- **Git**: For cloning the repository
- **C++ Compiler**: Required for optional PRIM analysis
- **RAM**: Minimum 4 GB recommended for full feature set

## Installation from Source

The recommended method for development and local use.

### Step 1: Clone the Repository

```bash
git clone https://github.com/gassraphael/AlphaPEM.git
cd AlphaPEM
```

### Step 2: Install Julia

**Linux / macOS:**
```bash
curl -fsSL https://install.julialang.org | sh
```

**Windows:**
```bash
winget install Julia -s msstore
```

### Step 3: Instantiate the Julia Environment

This downloads and precompiles all required dependencies.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

> **Note on HPC Clusters:**
    If you see a `GLFWError (API_UNAVAILABLE)` related to `GLMakie` on headless nodes, this is expected and safe to 
    ignore. AlphaPEM will automatically fall back to `CairoMakie` for rendering.

## Optional Features

### PRIM-Based Valid Parameter Region Analysis

Required only if you plan to use `examples/run_parameter_validity.jl` and the `Parametrisation.ValidParameterRegion`
module for identifying physically valid parameter ranges.

#### a) Install System Dependencies

**Linux (Debian/Ubuntu):**
```bash
sudo apt update && sudo apt install -y \
    r-base build-essential cmake \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libwebp-dev libpng-dev libtiff5-dev libjpeg-dev \
    libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
```

**macOS (Homebrew):**
```bash
brew install r
# Xcode Command Line Tools are automatically installed with full C++ support
```

**Windows:**
Download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) alongside R.

#### b) Clone the IRD Method Package

```bash
git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/IRD_method_2023
```

#### c) Install R Dependencies

```bash
sudo Rscript src/alphapem/parametrisation/validity/R/install_r_packages.R
```

On Windows, open Command Prompt as Administrator and omit `sudo`.

### Genetic Algorithm-Based Calibration

Required only for parameter calibration using `examples/run_calibration.jl`.

```bash
julia --project=. -e 'using CondaPkg; CondaPkg.resolve()'
```

This sets up Python/Conda dependencies automatically.

## Installation as a Package

To use AlphaPEM in other projects without cloning the full repository.

### From Julia

```julia
using Pkg
Pkg.add(url="https://github.com/gassraphael/AlphaPEM.git")
using AlphaPEM
```

### From Python

```bash
pip install juliacall
```

```python
from juliacall import Main as jl
jl.Pkg.add(url="https://github.com/gassraphael/AlphaPEM.git")
jl.seval("using AlphaPEM")
```

This approach allows seamless integration into polyglot applications with minimal overhead.

## Verify Installation

To verify that AlphaPEM is correctly installed and functional:

```bash
julia --project=. -e 'using AlphaPEM; println("AlphaPEM installed successfully!")'
```

## IDE Setup (Recommended)

### PyCharm

Install the [Flexible Julia](https://plugins.jetbrains.com/plugin/29356-flexible-julia) plugin for Julia support within PyCharm.

### VS Code

Install the official [Julia Language Support](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) extension.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `GLFWError` on cluster nodes | Expected behavior; AlphaPEM falls back to `CairoMakie` |
| Build fails with dependency errors | Run `julia --project=. -e 'using Pkg; Pkg.update()'` |
| IRD installation fails | Verify R installation with `Rscript --version` |
| Slow first run | First-time precompilation is expected; subsequent runs are fast |

For additional help, see the [Quick Start](quickstart.md) guide or [contact us](mailto:raphael.gass@univ-reunion.fr).
