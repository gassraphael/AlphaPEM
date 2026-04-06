# cell_node.jl
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
#   - `C_N2`   is carried by both GC nodes (always present in AnodeGCNode for
#               future flow-through-anode mode; always present in CathodeGCNode
#               for air operation).
#   - `eta_c`  (cathode overpotential) is localised at the CCL, where the
#               oxygen reduction reaction takes place.

# ────────────────────────────────────────────────────────────────────────────────
# Abstract root type
# ────────────────────────────────────────────────────────────────────────────────

abstract type AbstractCellNode end

# ────────────────────────────────────────────────────────────────────────────────
# Anode side
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at one anode gas-channel node."""
struct AnodeGCNode <: AbstractCellNode
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
end

"""Internal state at one anode gas-diffusion-layer node."""
struct AnodeGDLNode <: AbstractCellNode
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
end

"""Internal state at one anode microporous-layer node."""
struct AnodeMPLNode <: AbstractCellNode
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_H2 :: Float64   # Hydrogen concentration            (mol·m⁻³)
end

"""Internal state at the anode catalyst-layer node.
Contains `lambda` because the ionomer is present in this layer."""
struct AnodeCLNode <: AbstractCellNode
    T      :: Float64   # Temperature                     (K)
    C_v    :: Float64   # Water vapour concentration      (mol·m⁻³)
    s      :: Float64   # Liquid water saturation         (–)
    lambda :: Float64   # Ionomer water content           (–)
    C_H2   :: Float64   # Hydrogen concentration          (mol·m⁻³)
end

# ────────────────────────────────────────────────────────────────────────────────
# Electrolyte
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at the membrane node.
Only `T` and `lambda` are defined here: the membrane is impermeable to gases
and liquid water does not exist as a separate phase inside Nafion."""
struct MembraneNode <: AbstractCellNode
    T      :: Float64   # Temperature                     (K)
    lambda :: Float64   # Ionomer water content           (–)
end

# ────────────────────────────────────────────────────────────────────────────────
# Cathode side
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at the cathode catalyst-layer node.
Contains `lambda` (ionomer) and `eta_c` (overpotential of the ORR)."""
struct CathodeCLNode <: AbstractCellNode
    T      :: Float64   # Temperature                     (K)
    C_v    :: Float64   # Water vapour concentration      (mol·m⁻³)
    s      :: Float64   # Liquid water saturation         (–)
    lambda :: Float64   # Ionomer water content           (–)
    C_O2   :: Float64   # Oxygen concentration            (mol·m⁻³)
    eta_c  :: Float64   # Cathode overpotential           (V)
end

"""Internal state at one cathode microporous-layer node."""
struct CathodeMPLNode <: AbstractCellNode
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
end

"""Internal state at one cathode gas-diffusion-layer node."""
struct CathodeGDLNode <: AbstractCellNode
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
end

"""Internal state at one cathode gas-channel node."""
struct CathodeGCNode <: AbstractCellNode
    T    :: Float64   # Temperature                       (K)
    C_v  :: Float64   # Water vapour concentration        (mol·m⁻³)
    s    :: Float64   # Liquid water saturation           (–)
    C_O2 :: Float64   # Oxygen concentration              (mol·m⁻³)
    C_N2 :: Float64   # Nitrogen concentration            (mol·m⁻³)
end

# ────────────────────────────────────────────────────────────────────────────────
# Manifold nodes  (mixture state: pressure P, relative humidity Phi)
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state at one supply or exhaust manifold node.
Manifolds are characterised by mixture properties (P, Phi) rather than species concentrations.
"""
struct ManifoldNode <: AbstractCellNode
    P::Float64     # Pressure                            (Pa)
    Phi::Float64   # Relative humidity                   (–)
end

# ────────────────────────────────────────────────────────────────────────────────
# 1D MEA state  (one column = one GC node)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete 1D internal state for one gas-channel column.

Type parameters
---------------
nb_gdl : Int   Number of nodes in each GDL (anode and cathode share the same count).
nb_mpl : Int   Number of nodes in each MPL.

Fields follow the through-plane order from anode to cathode.
"""
struct MEAState1D{nb_gdl, nb_mpl}
    agc  :: AnodeGCNode
    agdl :: NTuple{nb_gdl, AnodeGDLNode}
    ampl :: NTuple{nb_mpl, AnodeMPLNode}
    acl  :: AnodeCLNode
    mem  :: MembraneNode
    ccl  :: CathodeCLNode
    cmpl :: NTuple{nb_mpl, CathodeMPLNode}
    cgdl :: NTuple{nb_gdl, CathodeGDLNode}
    cgc  :: CathodeGCNode
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
    nodes :: NTuple{nb_nodes, ManifoldNode}
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
# P2D fuel-cell state  (MEA stack only)
# ────────────────────────────────────────────────────────────────────────────────

"""Complete P2D (1D+1D) internal state of the fuel cell MEA stack.
Contains only the electrochemical core (gas channels along the flow direction).

Type parameters
---------------
nb_gdl : Int   Number of nodes per GDL.
nb_mpl : Int   Number of nodes per MPL.
nb_gc  : Int   Number of gas-channel nodes (spatial discretisation along the channel).
               Typical range: 1 – 10.

Fields:
  - `nodes`: MEA stack at each gas-channel position.
"""
struct FuelCellStateP2D{nb_gdl, nb_mpl, nb_gc}
    nodes :: NTuple{nb_gc, MEAState1D{nb_gdl, nb_mpl}}
end

