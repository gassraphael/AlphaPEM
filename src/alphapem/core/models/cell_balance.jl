# cell_balance.jl
#
# Typed representation of local conservation balances used in the ODE right-hand side.
# These structures are separate from state structures (see cell_node.jl):
# - state structs describe what the solver integrates;
# - balance structs describe temporary terms used to build d(state)/dt.

"""Local conservation terms for one transported quantity at one node.

The generic balance form is:
    accumulation ~ (J_in - J_out) / H + S

Units are intentionally context-dependent and must follow the existing physics
choices in the corresponding equation:
- J_in/J_out can be mol.m-2.s-1, kg.m-2.s-1, W.m-2, etc.
- S can be mol.m-3.s-1, kg.m-3.s-1, W.m-3, etc.
"""
struct NodeBalance
    J_in::Float64
    J_out::Float64
    S::Float64
end

NodeBalance() = NodeBalance(0.0, 0.0, 0.0)

abstract type AbstractLayerBalance end

# Anode side
struct AnodeGCBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    H2::NodeBalance
    N2::NodeBalance
    thermal::NodeBalance
end

struct AnodeGDLBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    H2::NodeBalance
    thermal::NodeBalance
end

struct AnodeMPLBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    H2::NodeBalance
    thermal::NodeBalance
end

struct AnodeCLBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    H2::NodeBalance
    lambda::NodeBalance
    thermal::NodeBalance
end

# Electrolyte
struct MembraneBalance <: AbstractLayerBalance
    lambda::NodeBalance
    thermal::NodeBalance
end

# Cathode side
# eta_c is intentionally excluded from NodeBalance; it follows a dedicated
# capacitive kinetics equation and will be handled separately.
struct CathodeCLBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    O2::NodeBalance
    lambda::NodeBalance
    thermal::NodeBalance
end

struct CathodeMPLBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    O2::NodeBalance
    thermal::NodeBalance
end

struct CathodeGDLBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    O2::NodeBalance
    thermal::NodeBalance
end

struct CathodeGCBalance <: AbstractLayerBalance
    vapour::NodeBalance
    liquid::NodeBalance
    O2::NodeBalance
    N2::NodeBalance
    thermal::NodeBalance
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold balances (mixture: P and Phi)
# ────────────────────────────────────────────────────────────────────────────────

"""Local conservation terms for pressure and relative humidity at a manifold node."""
struct ManifoldBalance <: AbstractLayerBalance
    P::NodeBalance      # Pressure balance
    Phi::NodeBalance    # Relative humidity balance
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D MEA balance  (one column = one GC node)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete 1D balance container for one gas-channel column."""
struct MEABalance1D{nb_gdl, nb_mpl}
    agc::AnodeGCBalance
    agdl::NTuple{nb_gdl, AnodeGDLBalance}
    ampl::NTuple{nb_mpl, AnodeMPLBalance}
    acl::AnodeCLBalance
    mem::MembraneBalance
    ccl::CathodeCLBalance
    cmpl::NTuple{nb_mpl, CathodeMPLBalance}
    cgdl::NTuple{nb_gdl, CathodeGDLBalance}
    cgc::CathodeGCBalance
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold line balances
# ────────────────────────────────────────────────────────────────────────────────

"""Complete balance container for one manifold line."""
struct ManifoldLineBalance{nb_nodes}
    nodes::NTuple{nb_nodes, ManifoldBalance}
end

# ────────────────────────────────────────────────────────────────────────────────
# P2D fuel-cell balances  (MEA stack only)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete P2D balance container for the MEA stack.
Manifold balances are kept separate (to be grouped later).
"""
struct FuelCellBalanceP2D{nb_gdl, nb_mpl, nb_gc}
    nodes :: NTuple{nb_gc, MEABalance1D{nb_gdl, nb_mpl}}
end

