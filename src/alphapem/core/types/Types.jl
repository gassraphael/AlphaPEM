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
# Fuel Cell specific types
include("fuelcell/cell_state.jl")
include("fuelcell/cell_derivative.jl")
include("fuelcell/cell_intermediates.jl")
include("fuelcell/cell_balance.jl")
include("fuelcell/auxiliary.jl")
include("fuelcell/cell_flows.jl")
include("fuelcell/simulation_outputs.jl")

# Electrolyzer specific types (to be added)
# include("electrolyzer/...")

# ── Public API ────────────────────────────────────────────────────────────────

# --- cell_state.jl ---
export AbstractCellState
export AnodeGCState, AnodeGDLState, AnodeMPLState, AnodeCLState
export MembraneState
export CathodeCLState, CathodeMPLState, CathodeGDLState, CathodeGCState
export ManifoldState
export CellState1D, ManifoldLine
export _ManifoldStateBundle, _ManifoldDerivativeBundle
export FuelCellStateP2D

# --- cell_derivative.jl ---
export AbstractCellDerivative
export AnodeGCDerivative, AnodeGDLDerivative, AnodeMPLDerivative, AnodeCLDerivative
export MembraneDerivative
export CathodeCLDerivative, CathodeMPLDerivative, CathodeGDLDerivative, CathodeGCDerivative
export ManifoldDerivative
export CellDerivative1D, ManifoldLineDerivative
export FuelCellDerivativeP2D
export MEADissolvedWaterDerivative, MEALiquidWaterDerivative, MEAVaporDerivative
export MEAGasSpeciesDerivative, MEAVoltageDerivative, MEATemperatureDerivative
export GCGasDerivative, GCLiquidWaterDerivative, GCTemperatureDerivative
export assemble_mea_derivative_1D, assemble_gc_derivative_1D

# --- cell_intermediates.jl ---
export AbstractCellIntermediate
export AnodeGCIntermediates, CathodeGCIntermediates
export MEAThermalIntermediates, CellIntermediates1D
export ManifoldNodeIntermediates, ManifoldIntermediates
export AuxiliaryIntermediates
export FuelCellIntermediatesP2D
export MEAFlowsIntWorkspace, MEAHeatIntWorkspace,
       GCManifoldWorkspace

# --- cell_balance.jl ---
export NodeBalance
export AbstractLayerBalance
export AnodeGCBalance, AnodeGDLBalance, AnodeMPLBalance, AnodeCLBalance
export MembraneBalance
export CathodeCLBalance, CathodeMPLBalance, CathodeGDLBalance, CathodeGCBalance
export ManifoldBalance
export CellBalance1D, ManifoldLineBalance
export FuelCellBalanceP2D

# --- auxiliary.jl ---
export Auxiliary0DState, Auxiliary0DDerivative

# --- cell_flows.jl ---
# scalar flow / source types
export MEASorptionSources, MEADissolvedWaterFlux, MEAWaterProductionSources
export MEAGasReactionSources
export MEAReactionHeat, MEASorptionHeat, MEAProtonHeat
# dimension-dependent MEA flux types
export MEALiquidFluxes, MEAVaporFluxes, MEAHydrogenFluxes, MEAOxygenFluxes, MEANitrogenFluxes
export MEALiquidSources, MEAVaporSources
export MEAThermalFluxes, MEALiquidHeat, MEAElectricHeat
# MEA top-level containers
export MEAFlows1D, MEAHeatFlows1D
export MEAFlowsWorkspace, MEAHeatWorkspace
# GC / along-channel flow types
export GCVaporFlows, GCLiquidFlows, GCHydrogenFlows, GCOxygenFlows, GCNitrogenFlows
export GCMassFlows, DesiredInletFlows
export GCManifoldFlows1D

# --- simulation_outputs.jl ---
export SolverTrajectory, DerivedOutputs, SimulationOutputs, FourierOutputs

end  # module Types

