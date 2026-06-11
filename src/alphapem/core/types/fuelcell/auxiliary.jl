# auxiliary.jl
#
# Typed representation of balance-of-plant (BoP) auxiliary 0D systems.
# These components are external to the electrochemical core (MEA + manifolds)
# and are modeled as point systems (zero-dimensional).

# ────────────────────────────────────────────────────────────────────────────────
# Auxiliary 0D systems state
# ────────────────────────────────────────────────────────────────────────────────

"""Internal state of auxiliary 0D systems (compressor, humidifiers, back-pressure valves).

These are modeled as point components (zero-dimensional). Each represents
a physical device in the balance-of-plant:

Fields:
  - Wcp: compressor mass flow rate                    (kg.s-1)
  - Wa_inj: anode humidifier water injection          (kg.s-1)
  - Wc_inj: cathode humidifier water injection        (kg.s-1)
  - Abp_a: anode back-pressure valve throttle area    (m²)
  - Abp_c: cathode back-pressure valve throttle area  (m²)
"""
struct Auxiliary0DState
    Wcp::Float64      # Compressor flow
    Wa_inj::Float64   # Anode humidifier injection
    Wc_inj::Float64   # Cathode humidifier injection
    Abp_a::Float64    # Anode back-pressure valve area
    Abp_c::Float64    # Cathode back-pressure valve area
end

# ────────────────────────────────────────────────────────────────────────────────
# Auxiliary 0D systems derivatives
# ────────────────────────────────────────────────────────────────────────────────

"""Time derivatives of auxiliary 0D systems.

Mirrors the structure of Auxiliary0DState for consistency with the ODE solver interface.
"""
struct Auxiliary0DDerivative
    Wcp::Float64
    Wa_inj::Float64
    Wc_inj::Float64
    Abp_a::Float64
    Abp_c::Float64
end
