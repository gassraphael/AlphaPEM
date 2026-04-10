# simulation_outputs.jl
#
# Typed simulation output containers used to replace string-keyed dictionaries
# in the runtime orchestration path.

"""Time trajectory of full typed solver states.

Type parameters
---------------
nb_gdl : Int   Number of GDL nodes per side.
nb_mpl : Int   Number of MPL nodes per side.
nb_gc  : Int   Number of GC nodes along the channel.
"""
struct SolverTrajectory{nb_gdl, nb_mpl, nb_gc}
    t      :: Vector{Float64}
    states :: Vector{FuelCellStateP2D{nb_gdl, nb_mpl, nb_gc}}
end

"""Derived simulation outputs used by post-processing in no-display workflows.

The vectors are stored per GC node when relevant.
"""
struct DerivedOutputs{nb_gc}
    Ucell  :: Vector{Float64}
    i_fc   :: Vector{Vector{Float64}}
    C_O2_Pt:: Vector{Vector{Float64}}
    v_a    :: Vector{Vector{Float64}}
    v_c    :: Vector{Vector{Float64}}
    Pa_in  :: Vector{Float64}
    Pc_in  :: Vector{Float64}
end

"""Complete typed outputs for one simulation run.

Groups solver trajectory and derived quantities while preserving compile-time
mesh information.
"""
struct SimulationOutputs{nb_gdl, nb_mpl, nb_gc}
    solver  :: SolverTrajectory{nb_gdl, nb_mpl, nb_gc}
    derived :: DerivedOutputs{nb_gc}
end

