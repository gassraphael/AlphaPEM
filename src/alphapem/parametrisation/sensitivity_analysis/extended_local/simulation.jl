# -*- coding: utf-8 -*-

# ============================================================
# SIMULATION RESULT
# ============================================================

"""
    ExtendedLocalSimulationResult

Results produced by local-sensitivity simulations, ready for downstream summarization and analysis.

# Fields
- `config::ExtendedLocalSensitivityConfig`: Configuration used for the runs.
- `params::Vector{ParameterBound}`: Analysed parameters in sensitivity-matrix row order.
- `screening_records::Vector{ScreeningRecord}`: Outcomes for all screened candidates.
- `candidates::Vector{ScreeningRecord}`: Accepted physical parameter candidates, in input order.
- `theta::Matrix{Float64}`: Accepted candidates mapped to the unit hypercube (`n_params × n_candidates`).
- `operating_cases::Vector{OperatingCase}`: Operating cases in simulation order.
- `sensitivities::Matrix{Union{Nothing,Matrix{Float64}}}`: One `n_params × n_current_densities`
  signed sensitivity block per `(candidate, operating case)`; incomplete blocks are `nothing`.
- `validity_masks::Matrix{Union{Nothing,BitMatrix}}`: Finite-sensitivity masks corresponding to each block.
- `statuses::Matrix{Symbol}`: Block status: `:complete`, `:partial`, or `:invalid`.
- `unit_times::Matrix{Float64}`: Wall-clock execution time per `(candidate, operating case)` (s).
- `run_status_counts::Dict{Symbol,Int}`: Counts of complete, partial, invalid, and thrown simulation runs.
- `run_errors::Dict{String,Int}`: Counts of caught simulation exceptions, grouped by error text.
- `execution_time::Float64`: Total local-analysis time, excluding warm-up (s).
"""
struct ExtendedLocalSimulationResult
    config::ExtendedLocalSensitivityConfig
    params::Vector{ParameterBound}
    screening_records::Vector{ScreeningRecord}
    candidates::Vector{ScreeningRecord}
    theta::Matrix{Float64}
    operating_cases::Vector{OperatingCase}
    sensitivities::Matrix{Union{Nothing, Matrix{Float64}}}
    validity_masks::Matrix{Union{Nothing, BitMatrix}}
    statuses::Matrix{Symbol}
    unit_times::Matrix{Float64}
    run_status_counts::Dict{Symbol, Int}
    run_errors::Dict{String, Int}
    execution_time::Float64
end


# ============================================================
# ALPHAPEM SIMULATOR AND WARM-UP
# ============================================================

"""
    alphapem_local_sensitivity_simulator(params, cfg) -> Function

Build the AlphaPEM simulation callable used by the local-sensitivity engine.

The callable accepts a unit-hypercube parameter vector and one operating case, runs AlphaPEM with
the configured ascending current profile, and returns the measured voltage curve.

Parameters
----------
params : Vector{ParameterBound}
    Analysed parameters in the same order as the unit-hypercube vector.
cfg : ExtendedLocalSensitivityConfig
    Fuel-cell, solver, and current-profile configuration.

Returns
-------
Function
    Callable `(theta, operating_case) -> VoltageCurve`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function alphapem_local_sensitivity_simulator(
    params::Vector{ParameterBound},
    cfg::ExtendedLocalSensitivityConfig,
)::Function
    # === Prepare values shared by every simulation ===
    # Validate once here because the returned closure may be called many times by the engine.
    validate_config(cfg)
    # Keep one nominal physical-parameter set so parameters not being analysed remain unchanged.
    base_fuel_cell = create_fuelcell(
        cfg.fuel_cell_type, cfg.voltage_zone; year = cfg.year, nb_gc = cfg.nb_gc)
    base_physical_parameters = base_fuel_cell.physical_parameters
    base_parameter_fields = (; (field => getfield(base_physical_parameters, field)
                                for field in fieldnames(PhysicalParams))...)
    # Reuse the solver settings and ascending current profile throughout this run.
    numerical_params = _numerical_params_for_local_simulation(cfg)
    current_params = ascending_polarization_params(cfg)

    return function (theta::AbstractVector{<:Real}, operating_case::OperatingCase)
        # === Validate and map one candidate's parameter coordinates ===
        length(theta) == length(params) ||
            throw(ArgumentError("theta must have one coordinate per analysed parameter"))

        # Convert normalized values to physical units, applying the configured integer policy.
        physical_values = [from_unit(
            parameter, theta[index]; round_int = cfg.integer_policy == :round)
                           for (index, parameter) in enumerate(params)]
        # Build field overrides using the types expected by PhysicalParams.
        overrides = (; (parameter.name =>
                        (parameter.type == :int ? Int(round(value)) : Float64(value))
                        for (parameter, value) in zip(params, physical_values))...)
        # Merge the analysed values with nominal values for every unanalysed physical parameter.
        physical_parameters = PhysicalParams(; merge(base_parameter_fields, overrides)...)

        # === Configure and run this candidate/case simulation ===
        # Operating conditions come from the case; the same current grid is used for all cases.
        simulation_config = SimulationConfig(
            type_fuel_cell = cfg.fuel_cell_type,
            year = cfg.year,
            type_current = current_params,
            numerical_parameters = numerical_params,
            voltage_zone = cfg.voltage_zone,
            type_display = :no_display,
            physical_parameters = physical_parameters,
            operating_conditions = operating_case.operating_conditions,
        )
        # run_simulation returns the simulator with solver outputs available for extraction.
        simulator = run_simulation(simulation_config)
        # Extract the measured voltage curve and its valid-current mask.
        return extract_voltage_curve(simulator)
    end
end


"""
    _numerical_params_for_local_simulation(cfg) -> NumericalParams

Copy the configured numerical mesh and tolerances while applying the local-run timeout.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Analysis configuration with mesh settings and maximum run time.

Returns
-------
NumericalParams
    Numerical settings passed to each AlphaPEM simulation.
"""
function _numerical_params_for_local_simulation(cfg::ExtendedLocalSensitivityConfig)::NumericalParams
    # Copy the configured mesh and tolerances; only the gas-channel count and run timeout
    # are taken from the extended-analysis configuration.
    np = cfg.numerical_params
    return NumericalParams(
        nb_gc = cfg.nb_gc,
        nb_gdl = np.nb_gdl,
        nb_mpl = np.nb_mpl,
        nb_man = np.nb_man,
        purge_time = np.purge_time,
        delta_purge = np.delta_purge,
        rtol = np.rtol,
        atol = np.atol,
        maxiters = np.maxiters,
        max_run_time_s = cfg.max_run_time_s,
        save_freq = np.save_freq,
    )
end


"""
    _warm_up_local_simulation(simulate_fn, n_params, operating_case)

Run and discard a midpoint simulation to compile the local simulation path before timing the run.

Parameters
----------
simulate_fn : Function
    Simulator callable for one unit-space parameter vector and operating case.
n_params : Int
    Number of analysed parameters, and therefore the dimension of the midpoint vector.
operating_case : OperatingCase
    First operating case used for compilation.

Returns
-------
Any
    The simulator's output, which the caller intentionally discards.
"""
function _warm_up_local_simulation(
    simulate_fn::Function,
    n_params::Int,
    operating_case::OperatingCase,
)
    # Compile the same call path used by real work units using a representative midpoint input.
    # Its output is discarded so it cannot be mistaken for an analysed candidate result.
    return simulate_fn(fill(0.5, n_params), operating_case)
end


# ============================================================
# ONE LOCAL-SENSITIVITY WORK UNIT
# ============================================================

"""
    _validate_local_voltage_curve(curve, current_density_grid) -> BitVector

Validate the extracted curve and return the finite, valid prefix mask.

Parameters
----------
curve : Any
    Value returned by the injected or AlphaPEM local simulator.
current_density_grid : Vector{Float64}
    Configured current grid in A/cm².

Returns
-------
BitVector
    Mask of finite valid points.

Raises
------
ArgumentError
    If the simulator returns a non-curve, a different grid, or a non-prefix mask.
"""
function _validate_local_voltage_curve(
    curve,
    current_density_grid::Vector{Float64},
)::BitVector
    # === Check the simulator's output contract ===
    # The current grid must match exactly so voltage columns have the same meaning in every run.
    curve isa VoltageCurve ||
        throw(ArgumentError("The local simulator must return a VoltageCurve"))
    length(curve.current_density) == length(current_density_grid) ||
        throw(ArgumentError("The local simulator returned the wrong number of current-density points"))
    all(isapprox.(curve.current_density, current_density_grid; atol = 1e-12, rtol = 0.0)) ||
        throw(ArgumentError("The local simulator returned a different current-density grid"))

    # Require both the profile's validity flag and a finite voltage for a usable measurement.
    valid = BitVector(curve.valid .& isfinite.(curve.voltage))
    # The ascending sweep cannot skip a failed current and later resume at a higher one.
    missing_seen = false
    for point_is_valid in valid
        missing_seen |= !point_is_valid
        point_is_valid && missing_seen &&
            throw(ArgumentError("A local voltage curve must have a valid-prefix mask"))
    end
    return valid
end


"""
    _run_local_sensitivity_unit(theta, operating_case, params, cfg, simulate_fn)

Compute all parameter sensitivities for one candidate and one operating case.

Each parameter is swept independently over the configured equidistant unit levels. A current point
contributes only when its voltage was measured in every included level used by the finite-difference
calculation. Missing points stay `NaN`; the sweep is never interpolated or extended past its valid
prefix.

Parameters
----------
theta : Vector{Float64}
    Accepted parameter candidate in unit-hypercube coordinates.
operating_case : OperatingCase
    Operating conditions for this simulation work unit.
params : Vector{ParameterBound}
    Analysed parameters in sensitivity row order.
cfg : ExtendedLocalSensitivityConfig
    Current grid, number of levels, and solver settings.
simulate_fn : Function
    Callable returning a `VoltageCurve` for one parameter vector and operating case.

Returns
-------
Tuple
    Sensitivity block, validity mask, block status, per-run status counts, and caught run errors.
"""
function _run_local_sensitivity_unit(
    theta::Vector{Float64},
    operating_case::OperatingCase,
    params::Vector{ParameterBound},
    cfg::ExtendedLocalSensitivityConfig,
    simulate_fn::Function,
)
    # === Initialize one sensitivity block and run-status diagnostics ===
    n_params = length(params)
    n_currents = length(cfg.current_density_grid)
    levels = equidistant_levels(cfg.n_levels)
    sensitivities = fill(NaN, n_params, n_currents)
    run_status_counts = Dict{Symbol, Int}(
        :complete => 0, :partial => 0, :invalid => 0, :failed => 0)
    run_errors = Dict{String, Int}()

    # === Evaluate parameter sweeps at this candidate and operating case ===
    # Sweep one parameter at a time; all other coordinates stay at the candidate values.
    for parameter_index in 1:n_params
        # Rows are parameter levels and columns are current densities. NaN marks a missing
        # or invalid voltage until a successful simulation fills that cell.
        voltage_by_level = fill(NaN, cfg.n_levels, n_currents)
        for (level_index, level) in enumerate(levels)
            # Copy the candidate so each simulation perturbs exactly one parameter.
            perturbed_theta = copy(theta)
            perturbed_theta[parameter_index] = level

            # Isolate a failed run and retain its message; cancellation and memory exhaustion
            # are process-level conditions and must not be converted into ordinary run failures.
            curve = try
                simulate_fn(perturbed_theta, operating_case)
            catch err
                err isa InterruptException && rethrow()
                err isa OutOfMemoryError && rethrow()
                error_key = string(nameof(typeof(err)), ": ", sprint(showerror, err))
                run_errors[error_key] = get(run_errors, error_key, 0) + 1
                run_status_counts[:failed] += 1
                continue
            end

            if curve === nothing
                # An explicit no-result return is an invalid run, but not a thrown exception.
                run_status_counts[:invalid] += 1
                continue
            end

            # Validate the return shape and grid separately so contract violations are diagnosable.
            valid = try
                _validate_local_voltage_curve(curve, cfg.current_density_grid)
            catch err
                error_key = string(nameof(typeof(err)), ": ", sprint(showerror, err))
                run_errors[error_key] = get(run_errors, error_key, 0) + 1
                run_status_counts[:failed] += 1
                continue
            end
            # Count a completed run separately from a useful partial prefix or an empty result.
            n_valid = count(valid)
            run_status_counts[n_valid == 0 ? :invalid :
                              (n_valid == n_currents ? :complete : :partial)] += 1
            for current_index in 1:n_currents
                valid[current_index] || continue
                voltage_by_level[level_index, current_index] = curve.voltage[current_index]
            end
        end

        # === Compute normalized sensitivities at each current density ===
        # A current contributes only if at least two parameter levels produced finite voltages.
        for current_index in 1:n_currents
            valid_levels = findall(isfinite, @view voltage_by_level[:, current_index])
            length(valid_levels) >= 2 || continue
            valid_values = voltage_by_level[valid_levels, current_index]
            # The output-mean normalization makes the sensitivity dimensionless.
            mean_voltage = mean(valid_values)
            isfinite(mean_voltage) && abs(mean_voltage) > eps(Float64) || continue
            # Use neighboring successful levels; when some levels failed, the level spacing
            # in the denominator reflects the wider interval between the remaining samples.
            slopes = diff(valid_values) ./ diff(levels[valid_levels])
            normalized_sensitivity = mean(slopes) / mean_voltage
            isfinite(normalized_sensitivity) &&
                (sensitivities[parameter_index, current_index] = normalized_sensitivity)
        end
    end

    # The block is complete only when every parameter/current entry is finite.
    # If none are finite it is invalid; otherwise it is partial.
    n_valid_sensitivities = count(isfinite, sensitivities)
    block_status = n_valid_sensitivities == 0 ? :invalid :
                   (n_valid_sensitivities == length(sensitivities) ? :complete : :partial)
    return sensitivities, BitMatrix(isfinite.(sensitivities)), block_status,
           run_status_counts, run_errors
end


# ============================================================
# CHECKPOINT IDENTITY AND PERSISTENCE
# ============================================================

"""
    _local_simulation_config_hash(cfg, params, screening_records, operating_cases, simulate_fn)

Hash all configuration and design values that can change local simulation results.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Configuration whose simulation-affecting fields are hashed.
params : Vector{ParameterBound}
    Analysed parameter definitions.
screening_records : AbstractVector{ScreeningRecord}
    Screening results, including accepted, rejected, and failed candidates.
operating_cases : Vector{OperatingCase}
    Ordered operating-case library.
simulate_fn : Function
    Simulator callable type included to distinguish injected simulation paths.

Returns
-------
String
    Lowercase SHA-256 digest of the local simulation design.
"""
function _local_simulation_config_hash(
    cfg::ExtendedLocalSensitivityConfig,
    params::Vector{ParameterBound},
    screening_records::AbstractVector{ScreeningRecord},
    operating_cases::Vector{OperatingCase},
    simulate_fn::Function,
)::String
    # Runtime/output controls do not change the mathematical result and are excluded by design.
    excluded_fields = Set((:output_dir, :parallel, :top_k, :resume, :checkpoint_frequency))
    parts = String[]
    # Serialize every remaining configuration field so changing a simulation setting invalidates
    # an old checkpoint instead of silently reusing results from a different run.
    for field in fieldnames(ExtendedLocalSensitivityConfig)
        field in excluded_fields && continue
        push!(parts, string(field, "=", _struct_to_string(getfield(cfg, field))))
    end
    append!(parts, ["parameter=" * _struct_to_string(parameter) for parameter in params])
    # Include rejected/failed records as well as accepted ones: they describe the screening result
    # set and preserve the identity and ordering of the candidate design.
    append!(parts, ["screening=$(record.id):$(record.status):$(_hash_sample(record.θ))"
                    for record in screening_records])
    # Operating conditions and their order determine each simulation's physical inputs.
    append!(parts, ["case=$(case.label):$(_struct_to_string(case.operating_conditions))"
                    for case in operating_cases])
    # Distinguish the default AlphaPEM simulator from an injected test/custom simulator.
    push!(parts, "simulator_type=$(typeof(simulate_fn))")
    return bytes2hex(SHA.sha256(codeunits(join(parts, "\n"))))
end


"""
    _save_local_sensitivity_checkpoint(path, config_hash, ...)

Atomically save completed sensitivity blocks and their design, timing, and diagnostic metadata.

Parameters
----------
path : String
    Destination JLD2 checkpoint path.
config_hash : String
    Digest of all simulation-affecting inputs.
candidate_ids, screening_ids, screening_statuses, screening_records
    Candidate identifiers and screening outcomes.
parameter_names, case_labels
    Ordered parameter and operating-case identifiers.
theta, sensitivities, validity_masks, statuses, unit_times, done
    Unit-space candidate matrix and completed block state.
run_status_counts, run_errors
    Aggregated simulation diagnostics.
elapsed_time : Float64
    Total elapsed analysis time excluding warm-up.

Returns
-------
Nothing
    The shared atomic writer reports persistence failures through the logger.
"""
function _save_local_sensitivity_checkpoint(
    path::String,
    config_hash::String,
    candidate_ids::Vector{Int},
    screening_ids::Vector{Int},
    screening_statuses::Vector{Symbol},
    screening_records::Vector{ScreeningRecord},
    parameter_names::Vector{Symbol},
    case_labels::Vector{String},
    theta::Matrix{Float64},
    sensitivities::Matrix{Union{Nothing, Matrix{Float64}}},
    validity_masks::Matrix{Union{Nothing, BitMatrix}},
    statuses::Matrix{Symbol},
    unit_times::Matrix{Float64},
    done::BitMatrix,
    run_status_counts::Dict{Symbol, Int},
    run_errors::Dict{String, Int},
    elapsed_time::Float64,
)
    _save_jld2_atomic(
        path,
        "extended local sensitivity";
        config_hash,
        candidate_ids,
        screening_ids,
        screening_statuses,
        screening_records,
        parameter_names,
        case_labels,
        theta,
        sensitivities = copy(sensitivities),
        validity_masks = copy(validity_masks),
        statuses = copy(statuses),
        unit_times = copy(unit_times),
        done = copy(done),
        run_status_counts = copy(run_status_counts),
        run_errors = copy(run_errors),
        elapsed_time,
    )
    return nothing
end


# ============================================================
# LOCAL-SENSITIVITY SIMULATION ENGINE
# ============================================================

"""
    run_local_sensitivity_simulations(records, params, operating_cases, cfg; simulate_fn=nothing)
        -> ExtendedLocalSimulationResult

Run the local sensitivity analysis for every accepted parameter candidate and operating case.

One work unit is `(accepted candidate, operating case)`. It sequentially runs `n_params × n_levels`
simulations, computes an `n_params × n_current_densities` sensitivity block, and retains a validity
mask. Independent work units are threaded when requested. A failed simulation is isolated and
reported in the returned diagnostics; partial current sweeps contribute only their valid prefix.

Parameters
----------
records : AbstractVector{ScreeningRecord}
    Screening outcomes; only records with status `:accepted` are simulated.
params : Vector{ParameterBound}
    Analysed parameters, in the same order as each candidate's physical parameter vector.
operating_cases : Vector{OperatingCase}
    Complete operating-case library.
cfg : ExtendedLocalSensitivityConfig
    Analysis and execution configuration.
simulate_fn : Function or nothing, optional
    Test/custom simulator with signature `(theta_unit, operating_case) -> VoltageCurve`. If
    `nothing`, use AlphaPEM with the configured physical parameters and current profile.

Returns
-------
ExtendedLocalSimulationResult
    Sensitivity blocks, masks, block statuses, current-case order, run diagnostics, and timing.

Raises
------
ArgumentError
    If no candidates are accepted, inputs have incompatible dimensions, or a checkpoint is invalid.
InterruptException
    Propagated so an interrupted run leaves its latest checkpoint available for resumption.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function run_local_sensitivity_simulations(
    records::AbstractVector{ScreeningRecord},
    params::Vector{ParameterBound},
    operating_cases::Vector{OperatingCase},
    cfg::ExtendedLocalSensitivityConfig;
    simulate_fn::Union{Nothing, Function} = nothing,
)::ExtendedLocalSimulationResult
    # === Validate inputs and select accepted candidates ===
    validate_config(cfg)
    isempty(params) && throw(ArgumentError("params must not be empty"))
    isempty(operating_cases) && throw(ArgumentError("operating_cases must not be empty"))
    candidate_ids_all = [record.id for record in records]
    length(unique(candidate_ids_all)) == length(candidate_ids_all) ||
        throw(ArgumentError("screening records must have unique candidate ids"))
    candidates = ScreeningRecord[record for record in records if record.status == :accepted]
    isempty(candidates) && throw(ArgumentError("No accepted screening candidates to analyse"))
    length(Set(parameter.name for parameter in params)) == length(params) ||
        throw(ArgumentError("params must have unique parameter names"))
    length(Set(case.label for case in operating_cases)) == length(operating_cases) ||
        throw(ArgumentError("operating_cases must have unique labels"))

    # === Convert accepted physical candidates to normalized coordinates ===
    # Each candidate stores physical values; the finite-difference analysis varies coordinates
    # in [0, 1], using the scale attached to each parameter.
    theta = Matrix{Float64}(undef, length(params), length(candidates))
    for (candidate_index, candidate) in enumerate(candidates)
        length(candidate.θ) == length(params) ||
            throw(ArgumentError("Candidate $(candidate.id) has the wrong number of parameter values"))
        for (parameter_index, parameter) in enumerate(params)
            coordinate = Float64(to_unit(parameter, candidate.θ[parameter_index]))
            isfinite(coordinate) && -1e-10 <= coordinate <= 1.0 + 1e-10 ||
                throw(ArgumentError("Candidate $(candidate.id) is outside the bounds for $(parameter.name)"))
            theta[parameter_index, candidate_index] = clamp(coordinate, 0.0, 1.0)
        end
    end

    # === Prepare the simulation callable and checkpoint identity ===
    # Tests can inject a deterministic callable; normal runs use the AlphaPEM-backed simulator.
    simulator = simulate_fn === nothing ?
        alphapem_local_sensitivity_simulator(params, cfg) : simulate_fn
    config_hash = _local_simulation_config_hash(
        cfg, params, records, operating_cases, simulator)
    checkpoint_path = joinpath(cfg.output_dir, "extended_local_sensitivity_checkpoint.jld2")

    # === Allocate candidate/case-indexed result and checkpoint storage ===
    # Matrix row = accepted candidate; matrix column = operating case.
    n_candidates = length(candidates)
    n_cases = length(operating_cases)
    sensitivities = Matrix{Union{Nothing, Matrix{Float64}}}(undef, n_candidates, n_cases)
    validity_masks = Matrix{Union{Nothing, BitMatrix}}(undef, n_candidates, n_cases)
    fill!(sensitivities, nothing)
    fill!(validity_masks, nothing)
    statuses = fill(:pending, n_candidates, n_cases)
    unit_times = fill(NaN, n_candidates, n_cases)
    done = falses(n_candidates, n_cases)
    candidate_ids = [candidate.id for candidate in candidates]
    screening_ids = [record.id for record in records]
    screening_statuses = [record.status for record in records]
    screening_records = collect(records)
    parameter_names = [parameter.name for parameter in params]
    case_labels = [operating_case.label for operating_case in operating_cases]
    restored_elapsed = 0.0
    # Keep run diagnostics alongside result blocks so a resumed run reports cumulative totals.
    run_status_counts = Dict{Symbol, Int}(
        :complete => 0, :partial => 0, :invalid => 0, :failed => 0)
    run_errors = Dict{String, Int}()

    # === Restore compatible checkpoint state when requested ===
    # A hash mismatch is rejected by the loader; matching metadata is checked again before
    # assigning saved blocks to their candidate/case positions.
    if cfg.resume
        checkpoint = _load_jld2_checked(
            checkpoint_path, config_hash, "extended local sensitivity")
        if checkpoint !== nothing
            checkpoint["candidate_ids"] == candidate_ids &&
            checkpoint["screening_ids"] == screening_ids &&
            checkpoint["screening_statuses"] == screening_statuses &&
            checkpoint["parameter_names"] == parameter_names &&
                checkpoint["case_labels"] == case_labels &&
                size(checkpoint["done"]) == size(done) ||
                throw(ArgumentError("Checkpoint metadata does not match the current local-sensitivity design"))
            done .= checkpoint["done"]
            statuses .= checkpoint["statuses"]
            unit_times .= checkpoint["unit_times"]
            restored_elapsed = checkpoint["elapsed_time"]
            merge!(run_status_counts, checkpoint["run_status_counts"])
            merge!(run_errors, checkpoint["run_errors"])
            for candidate_index in 1:n_candidates, case_index in 1:n_cases
                done[candidate_index, case_index] || continue
                # Only completed blocks are restored; unfinished entries remain `nothing`.
                block = checkpoint["sensitivities"][candidate_index, case_index]
                mask = checkpoint["validity_masks"][candidate_index, case_index]
                block isa Matrix{Float64} && mask isa BitMatrix ||
                    throw(ArgumentError("Checkpoint contains an invalid completed sensitivity block"))
                sensitivities[candidate_index, case_index] = block
                validity_masks[candidate_index, case_index] = mask
            end
            @info "Resuming extended local sensitivity" completed_units = count(done) total_units = length(done)
        end
    end

    # Always persist an initial checkpoint before warm-up or local simulation begins.
    # This leaves a valid resume file even if the first work unit is interrupted.
    _save_local_sensitivity_checkpoint(
        checkpoint_path, config_hash, candidate_ids, screening_ids, screening_statuses,
        screening_records,
        parameter_names, case_labels, theta, sensitivities, validity_masks,
        statuses, unit_times, done, run_status_counts, run_errors, restored_elapsed)

    # === Set up thread-safe progress, checkpoint, and diagnostic state ===
    lock_obj = ReentrantLock()
    warning_counters = _new_warning_counters()
    completed_since_checkpoint = 0
    pending_count = count(!, done)
    progress = pending_count == 0 ? nothing :
        Progress(pending_count; desc = "Extended local sensitivity: ", barlen = 40, color = :cyan)

    # === Limit BLAS threading and warm up Julia before timing ===
    # Outer work units provide parallelism, so one BLAS thread avoids nested oversubscription.
    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    execution_start = 0.0
    try
        return with_logger(_SimulationWarningLogger(warning_counters, current_logger())) do
            if simulate_fn === nothing && pending_count > 0
                @info "Pre-warming AlphaPEM local-sensitivity simulation..."
                try
                    _warm_up_local_simulation(simulator, length(params), first(operating_cases))
                catch err
                    err isa InterruptException && rethrow()
                    err isa OutOfMemoryError && rethrow()
                    @warn "AlphaPEM local-sensitivity warm-up failed; continuing with measured work units." exception = (err, catch_backtrace())
                end
            end

            # Warm-up time is excluded from analysis timing and the progress bar.
            execution_start = time()

            # Called only while holding `lock_obj`: snapshot arrays and cumulative diagnostics
            # together so the checkpoint cannot mix states from different completed units.
            function save_checkpoint_locked!()
                _save_local_sensitivity_checkpoint(
                    checkpoint_path, config_hash, candidate_ids, screening_ids,
                    screening_statuses, screening_records, parameter_names, case_labels, theta,
                    sensitivities, validity_masks, statuses, unit_times, done,
                    run_status_counts, run_errors,
                    restored_elapsed + time() - execution_start)
            end

            function run_one_unit!(candidate_index::Int, case_index::Int)
                # Keep the expensive simulations outside the lock; synchronize result commits,
                # shared diagnostics, checkpoints, and progress updates afterward.
                unit_start = time()
                block, mask, block_status, local_status_counts, local_errors =
                    _run_local_sensitivity_unit(
                        theta[:, candidate_index], operating_cases[case_index], params, cfg, simulator)
                elapsed = time() - unit_start

                # Commit this unit atomically with respect to checkpoint and progress updates.
                lock(lock_obj) do
                    sensitivities[candidate_index, case_index] = block
                    validity_masks[candidate_index, case_index] = mask
                    statuses[candidate_index, case_index] = block_status
                    unit_times[candidate_index, case_index] = elapsed
                    done[candidate_index, case_index] = true
                    for (status, count) in local_status_counts
                        run_status_counts[status] += count
                    end
                    for (message, count) in local_errors
                        run_errors[message] = get(run_errors, message, 0) + count
                    end

                    completed_since_checkpoint += 1
                    if cfg.checkpoint_frequency > 0 &&
                       completed_since_checkpoint % cfg.checkpoint_frequency == 0
                        save_checkpoint_locked!()
                    end
                    progress !== nothing && next!(progress)
                end
                return nothing
            end

            # === Execute independent candidate/case units serially or in parallel ===
            try
                if cfg.parallel && Threads.nthreads() > 1
                    # Spawn one task per unfinished matrix cell. Each task writes only its own cell;
                    # the lock protects shared counts, checkpoint snapshots, and progress updates.
                    @sync for candidate_index in 1:n_candidates, case_index in 1:n_cases
                        done[candidate_index, case_index] && continue
                        Threads.@spawn run_one_unit!(candidate_index, case_index)
                    end
                else
                    # Use the same unit function in serial mode to keep behavior identical.
                    for candidate_index in 1:n_candidates, case_index in 1:n_cases
                        done[candidate_index, case_index] && continue
                        run_one_unit!(candidate_index, case_index)
                    end
                end
            finally
                # Close the progress display even if a task throws or the run is interrupted.
                progress !== nothing && finish!(progress)
            end

            # === Report diagnostics and assemble the deterministic result ===
            elapsed = restored_elapsed + time() - execution_start
            @info "Extended local sensitivity units complete" complete_blocks = count(==(:complete), statuses) partial_blocks = count(==(:partial), statuses) invalid_blocks = count(==(:invalid), statuses) elapsed_seconds = round(elapsed; digits = 2)
            isempty(run_errors) ||
                @warn "Local simulation exceptions were isolated and recorded." error_counts = run_errors
            _print_warning_summary(warning_counters)

            # Assemble and return the completed simulation result in deterministic index order.
            result = ExtendedLocalSimulationResult(
                cfg, params, screening_records, candidates, theta, operating_cases,
                sensitivities, validity_masks, statuses, unit_times,
                run_status_counts, run_errors, elapsed)
            _remove_checkpoint(checkpoint_path)
            return result
        end
    finally
        # Restore the caller's BLAS setting even if warm-up, a task, or checkpointing fails.
        BLAS.set_num_threads(old_blas_threads)
    end
end


"""
    summarize_current_density_sensitivities(result, candidate_index, case_index) -> DataFrame

Summarize one completed sensitivity block from an `ExtendedLocalSimulationResult`.

Parameters
----------
result : ExtendedLocalSimulationResult
    Results returned by [`run_local_sensitivity_simulations`](@ref).
candidate_index : Int
    One-based position in `result.candidates`.
case_index : Int
    One-based position in `result.operating_cases`.

Returns
-------
DataFrame
    Full-grid and per-zone sensitivity magnitude, direction, and coverage summaries.

Raises
------
ArgumentError
    If the indices are out of range or the work unit has no completed block.
"""
function summarize_current_density_sensitivities(
    result::ExtendedLocalSimulationResult,
    candidate_index::Int,
    case_index::Int,
)::DataFrame
    # === Check the requested result indices and require a finished sensitivity block ===
    checkbounds(Bool, result.sensitivities, candidate_index, case_index) ||
        throw(ArgumentError("candidate_index or case_index is out of range"))
    block = result.sensitivities[candidate_index, case_index]
    block === nothing && throw(ArgumentError("The requested work unit has no completed sensitivity block"))
    # Delegate to the matrix-based summarizer so both entry points share identical zone rules.
    return summarize_current_density_sensitivities(block, result.params, result.config)
end
