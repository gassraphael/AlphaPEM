# cell_derivative.jl
#
# Typed representation of time derivatives used in the ODE right-hand side.
# These structures mirror the state topology from cell_state.jl.

abstract type AbstractCellDerivative end

# Anode side
struct AnodeGCDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    C_N2::Float64
    T::Float64
end

struct AnodeGDLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    T::Float64
end

struct AnodeMPLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    T::Float64
end

struct AnodeCLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_H2::Float64
    lambda::Float64
    T::Float64
end

# Electrolyte
struct MembraneDerivative <: AbstractCellDerivative
    lambda::Float64
    T::Float64
end

# Cathode side
struct CathodeCLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    lambda::Float64
    T::Float64
    eta_c::Float64
end

struct CathodeMPLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    T::Float64
end

struct CathodeGDLDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    T::Float64
end

struct CathodeGCDerivative <: AbstractCellDerivative
    C_v::Float64
    s::Float64
    C_O2::Float64
    C_N2::Float64
    T::Float64
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold derivatives
# ────────────────────────────────────────────────────────────────────────────────

"""Time derivatives for pressure and relative humidity at a manifold node."""
struct ManifoldDerivative <: AbstractCellDerivative
    P::Float64      # dP/dt
    Phi::Float64    # dPhi/dt
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D MEA derivative
# ────────────────────────────────────────────────────────────────────────────────

"""Complete 1D derivative container for one gas-channel column."""
struct MEACellDerivative1D{nb_gdl, nb_mpl}
    agc::AnodeGCDerivative
    agdl::NTuple{nb_gdl, AnodeGDLDerivative}
    ampl::NTuple{nb_mpl, AnodeMPLDerivative}
    acl::AnodeCLDerivative
    mem::MembraneDerivative
    ccl::CathodeCLDerivative
    cmpl::NTuple{nb_mpl, CathodeMPLDerivative}
    cgdl::NTuple{nb_gdl, CathodeGDLDerivative}
    cgc::CathodeGCDerivative
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold line derivatives
# ────────────────────────────────────────────────────────────────────────────────

"""Complete derivative container for one manifold line."""
struct ManifoldLineDerivative{nb_nodes}
    nodes::NTuple{nb_nodes, ManifoldDerivative}
end

# ────────────────────────────────────────────────────────────────────────────────
# P2D fuel-cell derivatives (MEA stack only)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete P2D derivative container for the MEA stack.
"""
struct FuelCellDerivativeP2D{nb_gdl, nb_mpl, nb_gc}
    nodes::NTuple{nb_gc, MEACellDerivative1D{nb_gdl, nb_mpl}}
end


