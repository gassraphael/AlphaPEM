# cell_state.jl
#
# Defines the Julia structs representing the internal state at each spatial node
# of the 1D MEA model and the 1D+1D (P2D) gas channel discretisation.
#
# Model layout (anode → cathode, along the through-plane direction):
#
#   AGC | AGDL×nb_gdl | AMPL×nb_mpl | ACL | Membrane | CCL | CMPL×nb_mpl | CGDL×nb_gdl | CGC
#
# Each struct holds only the state variables that are physically meaningful at
# that layer. In particular:
#   - `lambda` (ionomer water content) only appears in ACL, Membrane and CCL.
#   - `s`      (liquid water saturation) appears in all porous layers and GCs,
#               but the dominant transport mechanism differs:
#               convective in the GC, capillary-diffusive in GDL/MPL/CL.
#   - `C_N2`   is carried by both GC states (always present in AnodeGCState for
#               future flow-through-anode mode; always present in CathodeGCState
#               for air operation).
#   - `eta_c`  (cathode overpotential) is localised at the CCL, where the
#               oxygen reduction reaction takes place.

# ────────────────────────────────────────────────────────────────────────────────
# Abstract root type
# ────────────────────────────────────────────────────────────────────────────────

abstract type AbstractCellState end

# ────────────────────────────────────────────────────────────────────────────────
# Anode side
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at one anode gas-channel location."""
struct AnodeGCState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
end

"""Internal state at one anode gas-diffusion-layer location."""
struct AnodeGDLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
end

"""Internal state at one anode microporous-layer location."""
struct AnodeMPLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
end

"""Internal state at the anode catalyst-layer location.
Contains `lambda` because the ionomer is present in this layer."""
struct AnodeCLState <: AbstractCellState
    T      :: Float64   # Temperature                     (K)
    C_v    :: Float64   # Water vapour concentration      (mol·m⁻³)
    s      :: Float64   # Liquid water saturation         (–)
    lambda :: Float64   # Ionomer water content           (–)
    C_H2   :: Float64   # Hydrogen concentration          (mol·m⁻³)
end

# ────────────────────────────────────────────────────────────────────────────────
# Electrolyte
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at the membrane location.
Only `T` and `lambda` are defined here: the membrane is impermeable to gases
and liquid water does not exist as a separate phase inside Nafion."""
struct MembraneState <: AbstractCellState
    T      :: Float64   # Temperature                     (K)
    lambda :: Float64   # Ionomer water content           (–)
end

# ────────────────────────────────────────────────────────────────────────────────
# Cathode side
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at the cathode catalyst-layer location.
Contains `lambda` (ionomer) and `eta_c` (overpotential of the ORR)."""
struct CathodeCLState <: AbstractCellState
    T      :: Float64   # Temperature                     (K)
    C_v    :: Float64   # Water vapour concentration      (mol·m⁻³)
    s      :: Float64   # Liquid water saturation         (–)
    lambda :: Float64   # Ionomer water content           (–)
    C_O2   :: Float64   # Oxygen concentration            (mol·m⁻³)
    eta_c  :: Float64   # Cathode overpotential           (V)
end

"""Internal state at one cathode microporous-layer location."""
struct CathodeMPLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
end

"""Internal state at one cathode gas-diffusion-layer location."""
struct CathodeGDLState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
end

"""Internal state at one cathode gas-channel location."""
struct CathodeGCState <: AbstractCellState
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold nodes  (mixture state: pressure P, relative humidity Phi)
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at one supply or exhaust manifold location.
Manifolds are characterised by mixture properties (P, Phi) rather than species concentrations.
"""
struct ManifoldState <: AbstractCellState
    P::Float64     # Pressure                            (Pa)
    Phi::Float64   # Relative humidity                   (–)
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D cell-column state (MEA + AGC/CGC, one column = one GC node)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete 1D internal state for one gas-channel column.

Type parameters
---------------
nb_gdl : Int   Number of nodes in each GDL (anode and cathode share the same count).
nb_mpl : Int   Number of nodes in each MPL.

Fields follow the through-plane order from anode to cathode.
"""
struct CellState1D{nb_gdl, nb_mpl}
    agc  :: AnodeGCState
    agdl :: NTuple{nb_gdl, AnodeGDLState}
    ampl :: NTuple{nb_mpl, AnodeMPLState}
    acl  :: AnodeCLState
    mem  :: MembraneState
    ccl  :: CathodeCLState
    cmpl :: NTuple{nb_mpl, CathodeMPLState}
    cgdl :: NTuple{nb_gdl, CathodeGDLState}
    cgc  :: CathodeGCState
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold line (one per manifold: anode supply, anode exhaust, cathode supply, cathode exhaust)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete state of one manifold line (supply or exhaust).

Type parameters
---------------
nb_nodes : Int   Number of spatial nodes in this manifold.
                 Currently nb_nodes=1 (single-node manifold), but extensible.
"""
struct ManifoldLine{nb_nodes}
    nodes :: NTuple{nb_nodes, ManifoldState}
end

"""Typed bundle for the four manifold state lines.

This lightweight container keeps manifold entities grouped while preserving
separate manifold lines (asm, aem, csm, cem).
"""
struct _ManifoldStateBundle{ASM, AEM, CSM, CEM}
    asm::ASM
    aem::AEM
    csm::CSM
    cem::CEM
end

"""Typed bundle for the four manifold derivative lines.

Field types are generic to avoid coupling with declaration order between
node and equation model files.
"""
struct _ManifoldDerivativeBundle{ASM, AEM, CSM, CEM}
    asm::ASM
    aem::AEM
    csm::CSM
    cem::CEM
end
# ────────────────────────────────────────────────────────────────────────────────
# P2D fuel-cell state (cell columns = MEA + AGC/CGC)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete P2D (1D+1D) internal state of the fuel-cell cell-column stack.
Each node includes the full through-plane column (AGC + MEA core + CGC).

Type parameters
---------------
nb_gdl : Int   Number of nodes per GDL.
nb_mpl : Int   Number of nodes per MPL.
nb_gc  : Int   Number of gas-channel nodes (spatial discretisation along the channel).
               Typical range: 1 – 10.

Fields:
  - `nodes`: cell-column state at each gas-channel position.
"""
struct FuelCellStateP2D{nb_gdl, nb_mpl, nb_gc}
    nodes :: NTuple{nb_gc, CellState1D{nb_gdl, nb_mpl}}
end


