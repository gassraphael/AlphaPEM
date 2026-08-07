# -*- coding: utf-8 -*-

"""
    FuelCell factory for AlphaPEM

This module provides a factory function to instantiate a fuel cell from the
data stored under `data/`. No fuel-cell-specific source code is required:
adding a new stack only needs a new `data/<family>/<year>/stack.jl` file and
optional polarization YAML files under `data/<family>/<year>/pola/`.
"""

"""
    create_fuelcell(type_fuel_cell::Symbol, voltage_zone::Symbol; year::Union{Int,Nothing}=nothing)::AbstractFuelCell

Instantiate a fuel cell for `type_fuel_cell` from the data files.

The symbol is parsed to determine the `family` and `variant`. The factory then
loads `data/<family>/<year>/stack.jl` and, for variants,
`data/<family>/<year>/pola/<variant>.yml`. The `year` keyword is required when
the family is versioned by year (e.g. `data/ZSW/<year>/`); it is ignored for
families that are not versioned (e.g. `data/EH31/`).
"""
function create_fuelcell(type_fuel_cell::Symbol, voltage_zone::Symbol; year::Union{Int,Nothing}=nothing)::AbstractFuelCell
    family, variant = ExperimentalDataLoader.parse_fuel_cell_symbol(type_fuel_cell)

    # Provide a default year for families that are versioned by year.
    # Families that do not use a year directory ignore this value.
    if isnothing(year) && family == :ZSW
        year = 2024
    end

    return _load_fuelcell_from_data(family, variant, voltage_zone, year)
end

"""
    _load_fuelcell_from_data(family::Symbol, variant::String, voltage_zone::Symbol,
                             year::Union{Int,Nothing}=nothing)::AbstractFuelCell

Generic constructor: load the stack file and the requested polarization data,
then return a `GenericFuelCell` populated with the loaded values.
"""
function _load_fuelcell_from_data(family::Symbol, variant::String, voltage_zone::Symbol,
                                  year::Union{Int,Nothing}=nothing)::AbstractFuelCell
    stack = ExperimentalDataLoader.load_stack_data(family, year)

    # Operating conditions: nominal from stack.jl, variant-specific from pola YAML.
    if variant == "nominal"
        oc = stack.operating_conditions
    else
        oc = ExperimentalDataLoader.load_operating_conditions(family, variant, year)
    end

    # Apply stack-specific corrections (e.g. ZSW +3 °C temperature correction).
    # `invokelatest` is required because the correction function is defined in a
    # dynamically loaded module and would otherwise be too new for this world age.
    oc = Base.invokelatest(stack.apply_operating_conditions_correction, oc)

    # Polarization data: nominal also has a dedicated YAML file.
    pola_exp, pola_exp_cali = ExperimentalDataLoader.load_pola_experimental_data(
        family, variant, voltage_zone, year)

    undet = Dict{Symbol, Vector{Tuple{Symbol, Float64, Float64}}}(
        :full => stack.undetermined_parameters_full,
        :before_voltage_drop => stack.undetermined_parameters_before_voltage_drop,
    )

    return GenericFuelCell(
        physical_parameters = stack.physical_parameters,
        operating_conditions = oc,
        pola_exp_data = pola_exp,
        pola_exp_data_cali = pola_exp_cali,
        undetermined_parameters = undet
    )
end
