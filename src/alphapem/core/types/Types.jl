"""
    AlphaPEM.Core.Types

Pure data types (structs) for the PEMFC domain model.
These definitions are dependency-free (Base only) and must be loaded before
any computational module that constructs or dispatches on these types.

Structs are grouped by physical role:
  - cell_state        : integration variables at each spatial node
  - cell_derivative   : time derivatives (ODE RHS outputs)
  - cell_intermediates: scratch values computed during the ODE RHS evaluation
  - cell_balance      : local conservation-law terms (fluxes + sources)
  - auxiliary         : balance-of-plant 0-D auxiliary systems
  - cell_flows        : inter-layer fluxes and source terms (ODE RHS outputs)
"""
module Types

# ── Type definition files (no imports required – Base types only) ────────────
include("cell_state.jl")
include("cell_derivative.jl")
include("cell_intermediates.jl")
include("cell_balance.jl")
include("auxiliary.jl")
include("cell_flows.jl")

# ── Public API ────────────────────────────────────────────────────────────────

# --- cell_state.jl ---
export AbstractCellState
export AnodeGCState, AnodeGDLState, AnodeMPLState, AnodeCLState
export MembraneState
export CathodeCLState, CathodeMPLState, CathodeGDLState, CathodeGCState
export ManifoldState
export MEAState1D, ManifoldLine
export _ManifoldStateBundle, _ManifoldDerivativeBundle
export FuelCellStateP2D

# --- cell_derivative.jl ---
export AbstractCellDerivative
export AnodeGCDerivative, AnodeGDLDerivative, AnodeMPLDerivative, AnodeCLDerivative
export MembraneDerivative
export CathodeCLDerivative, CathodeMPLDerivative, CathodeGDLDerivative, CathodeGCDerivative
export ManifoldDerivative
export MEACellDerivative1D, ManifoldLineDerivative
export FuelCellDerivativeP2D
export MEADissolvedWaterDerivative, MEALiquidWaterDerivative, MEAVaporDerivative
export MEAH2O2SpeciesDerivative, MEAVoltageDerivative, MEATemperatureDerivative

# --- cell_intermediates.jl ---
export AbstractCellIntermediate
export AnodeGCIntermediates, CathodeGCIntermediates
export MEAThermalIntermediates, MEAIntermediates1D
export ManifoldNodeIntermediates, ManifoldIntermediates
export AuxiliaryIntermediates
export FuelCellIntermediatesP2D

# --- cell_balance.jl ---
export NodeBalance
export AbstractLayerBalance
export AnodeGCBalance, AnodeGDLBalance, AnodeMPLBalance, AnodeCLBalance
export MembraneBalance
export CathodeCLBalance, CathodeMPLBalance, CathodeGDLBalance, CathodeGCBalance
export ManifoldBalance
export MEABalance1D, ManifoldLineBalance
export FuelCellBalanceP2D

# --- auxiliary.jl ---
export Auxiliary0DState, Auxiliary0DDerivative

# --- cell_flows.jl ---
# scalar flow / source types
export MEASorptionSources, MEADissolvedWaterFlux, MEAWaterProductionSources
export MEAGasReactionSources
export MEAReactionHeat, MEASorptionHeat, MEAProtonHeat
# dimension-dependent flux types
export MEALiquidFluxes, MEAVaporFluxes, MEAHydrogenFluxes, MEAOxygenFluxes
export MEALiquidSources, MEAVaporSources
export MEAThermalFluxes, MEALiquidHeat, MEAElectricHeat
# top-level containers
export MEAFlows1D, MEAHeatFlows1D

end  # module Types

