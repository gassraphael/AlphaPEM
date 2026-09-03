# -*- coding: utf-8 -*-

"""
    FuelCell factory for AlphaPEM

This module provides a factory function to instantiate a fuel cell from the
data stored under `data/`. No fuel-cell-specific source code is required:
adding a new stack only needs a new `data/<family>/<year>/stack.jl` file and
optional polarization YAML files under `data/<family>/<year>/pola/`.
"""

"""
    create_fuelcell(type_fuel_cell::Symbol, voltage_zone::Symbol; year::Union{Int,Nothing}=nothing, nb_gc::Int=5)::AbstractFuelCell

Instantiate a fuel cell for `type_fuel_cell` from the data files.

The symbol is parsed to determine the `family` and `variant` from directory and
YAML file names (e.g. `:EH_nominal`, `:ZSW_Pa_1_61_Pc_1_41`). The factory then
loads `data/<family>/<year>/stack.jl` and, for variants,
`data/<family>/<year>/pola/<variant>.yml`. The `year` keyword is required when
the family is versioned by year (e.g. `data/ZSW/<year>/`); it is ignored for
families that are not versioned (e.g. `data/EH/`).

The `nb_gc` keyword selects the physical-parameter set and undetermined-parameter
bounds stored in the returned fuel cell (`:1D` for `nb_gc = 1`, `:1D1D` otherwise).
"""
function create_fuelcell(type_fuel_cell::Symbol, voltage_zone::Symbol; year::Union{Int,Nothing}=nothing, nb_gc::Int=5)::AbstractFuelCell
    family, variant = ExperimentalDataLoader.parse_fuel_cell_symbol(type_fuel_cell)

    if !isdir(joinpath(ExperimentalDataLoader.project_root(), "data", string(family)))
        throw(ArgumentError("Unsupported fuel_cell_type: $type_fuel_cell"))
    end

    # The year is required for families that are versioned by year.
    # Families that do not use a year directory ignore this value.
    if isnothing(year) && isdir(joinpath(ExperimentalDataLoader.project_root(), "data", string(family)))
        # Check whether the family has year subdirectories and pick the most
        # recent one that contains a stack.jl as a default. This keeps the API
        # generic for future versioned families without hard-coding a year.
        data_dir = joinpath(ExperimentalDataLoader.project_root(), "data", string(family))
        year_dirs = [d for d in readdir(data_dir)
                     if isdir(joinpath(data_dir, d)) && all(isdigit, d) &&
                        isfile(joinpath(data_dir, d, "stack.jl")) &&
                        isdir(joinpath(data_dir, d, "pola"))]
        if !isempty(year_dirs)
            year = parse(Int, maximum(year_dirs))
        end
    end

    return _load_fuelcell_from_data(family, variant, voltage_zone, year; nb_gc=nb_gc)
end

"""
    _load_fuelcell_from_data(family::Symbol, variant::String, voltage_zone::Symbol,
                             year::Union{Int,Nothing}=nothing; nb_gc::Int=5)::AbstractFuelCell

Generic constructor: load the stack file and the requested polarization data,
then return a `GenericFuelCell` populated with the loaded values for the
requested GC spatial resolution.
"""
function _load_fuelcell_from_data(family::Symbol, variant::String, voltage_zone::Symbol,
                                  year::Union{Int,Nothing}=nothing; nb_gc::Int=5)::AbstractFuelCell
    stack = ExperimentalDataLoader.load_stack_data(family, year)

    # Operating conditions: nominal from stack.jl, variant-specific from pola YAML
    # if the polarization file exists.
    path = ExperimentalDataLoader.pola_file(family, variant, year)
    if variant == "nominal" || !isfile(path)
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

    discr = nb_gc == 1 ? Symbol("1D") : Symbol("1D1D")
    pp = stack.physical_parameters[discr]
    undet = stack.undetermined_parameter_bounds[discr]

    return GenericFuelCell(
        physical_parameters = pp,
        operating_conditions = oc,
        pola_exp_data = pola_exp,
        pola_exp_data_cali = pola_exp_cali,
        undetermined_parameter_bounds = undet
    )
end
