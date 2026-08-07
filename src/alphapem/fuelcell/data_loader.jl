# src/alphapem/fuelcell/data_loader.jl

"""
    ExperimentalDataLoader

Utilities for loading fuel-cell experimental data from YAML files under `data/`.

The source files keep current densities in A/cm²; they are converted to A/m²
when building `PolaExperimentalData` to match the internal SI convention.
"""
module ExperimentalDataLoader

using YAML

using ...Config: PhysicalParams, OperatingConditions, PolaExperimentalData

export project_root,
       fuel_cell_data_dir,
       stack_file,
       pola_file,
       load_stack_data,
       load_pola_data,
       load_operating_conditions,
       load_pola_experimental_data,
       build_operating_conditions,
       build_pola_experimental_data,
       parse_fuel_cell_symbol

# ═══════════════════════════════════════════════════════════════════════════════
# Project root resolution
# ═══════════════════════════════════════════════════════════════════════════════

"""
    project_root() -> String

Return the AlphaPEM project root directory.

The root is located by walking up from this file until a `Project.toml` marker
is found. This guarantees a stable path regardless of the current working
directory.
"""
function project_root()::String
    markers = ["Project.toml", "Manifest.toml", ".git"]
    parent = dirname(abspath(@__FILE__))
    while true
        any(m -> ispath(joinpath(parent, m)), markers) && return parent
        new_parent = dirname(parent)
        new_parent == parent && return pwd()
        parent = new_parent
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Fuel-cell symbol parsing
# ═══════════════════════════════════════════════════════════════════════════════

"""
    parse_fuel_cell_symbol(type_fuel_cell::Symbol) -> (family::Symbol, variant::String)

Parse a fuel-cell symbol into a `(family, variant)` pair used to locate data
files.

The symbol is expected to follow the convention `<FAMILY>_<STACK>[_<VARIANT>]`.
The family is the first underscore-separated token and is used as the directory
name under `data/`. The variant is everything after the stack name; if absent,
it defaults to `"nominal"`.

# Examples
- `:ZSW_GenStack`              -> `(:ZSW, "nominal")`
- `:ZSW_GenStack_T_84`         -> `(:ZSW, "T_84")`
- `:ZSW_GenStack_Pa_2_8_Pc_2_6` -> `(:ZSW, "Pa_2_8_Pc_2_6")`
- `:EH31_2022` / `:EH31_1_5`  -> `(:EH31, "nominal")`
- `:EH31_2_0`                 -> `(:EH31, "Pa_2_0")`
"""
function parse_fuel_cell_symbol(type_fuel_cell::Symbol)::Tuple{Symbol, String}
    s = string(type_fuel_cell)

    # The family is the longest leading prefix that matches an existing data
    # directory under `data/`. This makes the parser generic: any new family
    # only needs a new directory, with no source-code change.
    family = _detect_family_prefix(s)

    family_str = string(family)
    rest = startswith(s, family_str * "_") ? s[length(family_str)+2:end] : s[length(family_str)+1:end]

    # The stack name is the next token (e.g. "GenStack" in "ZSW_GenStack_T_84").
    # Everything after the stack name is the variant.
    stack_tokens = split(rest, '_')
    if rest == "" || length(stack_tokens) <= 1
        variant = "nominal"
    else
        variant = join(stack_tokens[2:end], "_")
    end

    return family, variant
end

"""
    _detect_family_prefix(s::String) -> Symbol

Return the longest leading prefix of `s` that matches an existing directory
under `data/`. Falls back to the first token if no directory matches.
"""
function _detect_family_prefix(s::String)::Symbol
    tokens = split(s, '_')
    data_root = joinpath(project_root(), "data")

    # Try progressively longer prefixes.
    for n in length(tokens):-1:1
        candidate = join(tokens[1:n], '_')
        isdir(joinpath(data_root, candidate)) && return Symbol(candidate)
    end

    # Fallback to the first token.
    return Symbol(tokens[1])
end

# ═══════════════════════════════════════════════════════════════════════════════
# Path helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fuel_cell_data_dir(family::Symbol, year::Union{Int,Nothing}=nothing) -> String

Return the directory that holds the data for a given fuel-cell family.

If `year` is provided, the directory is `data/<family>/<year>/`. Otherwise it is
`data/<family>/`. The function checks that the resulting directory exists.
"""
function fuel_cell_data_dir(family::Symbol, year::Union{Int,Nothing}=nothing)::String
    root = project_root()
    dir = if isnothing(year)
        joinpath(root, "data", string(family))
    else
        joinpath(root, "data", string(family), string(year))
    end
    isdir(dir) || throw(ArgumentError("Fuel-cell data directory not found: $dir"))
    return dir
end

"""
    stack_file(family::Symbol, year::Union{Int,Nothing}=nothing) -> String

Return the path to the master stack file (`stack.jl`).
"""
function stack_file(family::Symbol, year::Union{Int,Nothing}=nothing)::String
    return joinpath(fuel_cell_data_dir(family, year), "stack.jl")
end

"""
    pola_file(family::Symbol, variant::String, year::Union{Int,Nothing}=nothing) -> String

Return the path to a polarization data YAML file.
"""
function pola_file(family::Symbol, variant::String, year::Union{Int,Nothing}=nothing)::String
    return joinpath(fuel_cell_data_dir(family, year), "pola", variant * ".yml")
end

# ═══════════════════════════════════════════════════════════════════════════════
# YAML loading and conversion
# ═══════════════════════════════════════════════════════════════════════════════

"""
    load_yaml(path::String) -> Dict

Load a YAML file and return it as a `Dict`.
"""
function load_yaml(path::String)::Dict{Any, Any}
    isfile(path) || throw(ArgumentError("Data file not found: $path"))
    return YAML.load_file(path)
end

"""
    build_operating_conditions(oc::Dict) -> OperatingConditions

Build an `OperatingConditions` instance from a dictionary of field values.

Temperatures in the YAML files are stored in degrees Celsius and converted to
Kelvin here to match the internal SI convention.
"""
function build_operating_conditions(oc::Dict)::OperatingConditions
    return OperatingConditions(
        T_des       = Float64(oc["T_des"]) + 273.15,
        Pa_des      = Float64(oc["Pa_des"]),
        Pc_des      = Float64(oc["Pc_des"]),
        Sa          = Float64(oc["Sa"]),
        Sc          = Float64(oc["Sc"]),
        Phi_a_des   = Float64(oc["Phi_a_des"]),
        Phi_c_des   = Float64(oc["Phi_c_des"]),
        y_H2_in     = Float64(oc["y_H2_in"]),
        i_min_stoich= Float64(oc["i_min_stoich"]),
    )
end

"""
    build_pola_experimental_data(pola::Dict) -> PolaExperimentalData

Build a `PolaExperimentalData` instance from a dictionary.

Current densities in the YAML are stored in A/cm² and are converted to A/m²
by multiplying by 1e4.
"""
function build_pola_experimental_data(pola::Dict)::PolaExperimentalData
    i_exp = Float64.(pola["i_exp"]) .* 1e4
    U_exp = Float64.(pola["U_exp"])
    return PolaExperimentalData(i_exp = i_exp, U_exp = U_exp)
end

"""
    load_stack_data(family::Symbol, year::Union{Int,Nothing}=nothing)

Load the master stack file (`stack.jl`) for the requested fuel-cell family.

The stack file must define:
  - `PHYSICAL_PARAMETERS::PhysicalParams`
  - `OPERATING_CONDITIONS::OperatingConditions`
  - `UNDETERMINED_PARAMETERS_FULL::Vector{Tuple{Symbol, Float64, Float64}}`
  - `UNDETERMINED_PARAMETERS_BEFORE_VOLTAGE_DROP::Vector{Tuple{Symbol, Float64, Float64}}`

Returns a `NamedTuple` with the same fields.
"""
function load_stack_data(family::Symbol, year::Union{Int,Nothing}=nothing)
    path = stack_file(family, year)
    mod = Module()
    # Make AlphaPEM available so that `using AlphaPEM.Config` works inside stack.jl
    Core.eval(mod, :(using AlphaPEM.Config: PhysicalParams, OperatingConditions))
    Base.include(mod, path)

    correction = if isdefined(mod, :apply_operating_conditions_correction)
        Core.eval(mod, :apply_operating_conditions_correction)
    else
        identity
    end

    return (
        physical_parameters       = Core.eval(mod, :PHYSICAL_PARAMETERS),
        operating_conditions      = Core.eval(mod, :OPERATING_CONDITIONS),
        undetermined_parameters_full = Core.eval(mod, :UNDETERMINED_PARAMETERS_FULL),
        undetermined_parameters_before_voltage_drop = Core.eval(mod, :UNDETERMINED_PARAMETERS_BEFORE_VOLTAGE_DROP),
        apply_operating_conditions_correction = correction,
    )
end

"""
    load_pola_data(family::Symbol, variant::String, year::Union{Int,Nothing}=nothing) -> Dict

Load a polarization YAML file and return its parsed content.
"""
function load_pola_data(family::Symbol, variant::String, year::Union{Int,Nothing}=nothing)::Dict
    return load_yaml(pola_file(family, variant, year))
end

"""
    load_pola_experimental_data(family::Symbol, variant::String, voltage_zone::Symbol,
                                year::Union{Int,Nothing}=nothing)
        -> Tuple{PolaExperimentalData, PolaExperimentalData}

Load the polarization YAML file for `variant` and return the full and calibration
experimental data for the requested `voltage_zone`.

`voltage_zone` must be `:full` or `:before_voltage_drop`.
"""
function load_pola_experimental_data(family::Symbol, variant::String, voltage_zone::Symbol,
                                     year::Union{Int,Nothing}=nothing)::Tuple{PolaExperimentalData, PolaExperimentalData}
    voltage_zone in (:full, :before_voltage_drop) ||
        throw(ArgumentError("voltage_zone must be :full or :before_voltage_drop (got $voltage_zone)"))

    data = load_pola_data(family, variant, year)
    zone_key = string(voltage_zone)

    pola_exp = build_pola_experimental_data(data["pola_exp_data"][zone_key])
    pola_exp_cali = build_pola_experimental_data(data["pola_exp_data_calibration"][zone_key])

    return pola_exp, pola_exp_cali
end

"""
    load_operating_conditions(family::Symbol, variant::String, year::Union{Int,Nothing}=nothing)
        -> OperatingConditions

Load the operating conditions stored in a polarization YAML file.
"""
function load_operating_conditions(family::Symbol, variant::String, year::Union{Int,Nothing}=nothing)::OperatingConditions
    data = load_pola_data(family, variant, year)
    return build_operating_conditions(data["operating_conditions"])
end

end  # module ExperimentalDataLoader
