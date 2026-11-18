# -*- coding: utf-8 -*-

"""
    generate_unit_design(n_candidates, d; seed, nominal_unit = nothing) -> Matrix{Float64}

Generate an Owen-scrambled Sobol design in the unit hypercube.

A Sobol sequence is a **quasi-random** (quasi-Monte Carlo) sequence: unlike independent
pseudo-random draws, its points are constructed to cover the unit hypercube evenly, with low
discrepancy. This is the type of design meant by Goshtasbi et al.'s phrase "quasi-random Sobol
sequence." Owen scrambling randomizes the Sobol points while retaining their low-discrepancy
structure. Thus, this implementation uses a randomized version of the same Sobol design, not a
different sampling method. The `seed` makes that randomization reproducible; the cited phrase alone
does not specify whether scrambling was used in the article.

If `n_candidates` is not a power of 2, it is automatically rounded up to the next power of 2 and an
informational message reports the effective count. This preserves the Sobol balance guarantee.
Scrambling is always applied and seeded by `seed`, so the same seed produces the same design.
If `nominal_unit` (a `d`-vector in `[0, 1]`) is given, it is prepended as the first column, so the
design then has one more column than the effective Sobol count.

Parameters
----------
n_candidates : Int
    Requested number of Sobol candidates; rounded up to a power of 2 when necessary.
d : Int
    Number of dimensions of the unit hypercube.
seed : Int
    Seed used for Owen scrambling.
nominal_unit : AbstractVector{<:Real} or nothing, optional
    Nominal point to prepend to the design. Default is `nothing`.

Returns
-------
Matrix{Float64}
    Design with one candidate per column; includes the nominal point as the first column when supplied.

Raises
------
ArgumentError
    If `n_candidates` or `d` is not positive, or `nominal_unit` has the wrong dimension.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function generate_unit_design(n_candidates::Int, d::Int;
                              seed::Int,
                              nominal_unit::Union{Nothing, AbstractVector{<:Real}} = nothing)::Matrix{Float64}
    # Preliminary checks
    n_candidates >= 1 || throw(ArgumentError("n_candidates must be positive"))
    d >= 1 || throw(ArgumentError("d must be positive"))

    # Rounds the number of candidates up to the next power of 2.
    n_effective = nextpow(2, n_candidates)
    if n_effective != n_candidates
        @info "Requested $n_candidates Sobol candidates; rounding up to $n_effective to preserve the balance of a power-of-two design."
    end

    # Generate the scrambled Sobol design
    rng = Random.MersenneTwister(seed) # seeded pseudo-random number generator used to randomize the Sobol sequence.
    sampler = QuasiMonteCarlo.SobolSample(QuasiMonteCarlo.OwenScramble(; base = 2, rng = rng)) # generate Sobol points with an Owen scramble driven by `rng`.
    design = Matrix{Float64}(QuasiMonteCarlo.sample(n_effective, d, sampler)) # sampled unit-cube points: one point per column and one parameter per row.

    # Optionally prepend the nominal point
    if nominal_unit !== nothing
        length(nominal_unit) == d || throw(ArgumentError("nominal_unit must have length $d"))
        design = hcat(Float64.(nominal_unit), design)
    end

    return design
end


"""
    nominal_unit_point(params, physical_parameters) -> Vector{Float64}

Map the stack's nominal physical parameters into the unit hypercube.

Parameters
----------
params : Vector{ParameterBound}
    Parameter bounds and scales defining the unit-hypercube coordinates.
physical_parameters : PhysicalParams
    Nominal physical parameters of the stack.

Returns
-------
Vector{Float64}
    Nominal point in the unit hypercube, clamped coordinate-wise to `[0, 1]`.
"""
function nominal_unit_point(params::Vector{ParameterBound}, physical_parameters)::Vector{Float64}
    return [clamp(Float64(to_unit(p, Float64(getfield(physical_parameters, p.name)))), 0.0, 1.0)
            for p in params]
end


"""
    unit_design_to_physical(design, params) -> Matrix{Float64}

Map a unit-hypercube design to physical parameter values.

Each row is transformed with its parameter's linear or logarithmic scale, and integer parameters are
rounded to the nearest integer.

Parameters
----------
design : AbstractMatrix{<:Real}
    Unit-hypercube design with one candidate per column.
params : Vector{ParameterBound}
    Bounds and scales for the design dimensions.

Returns
-------
Matrix{Float64}
    Physical design with the same shape as `design`.

Raises
------
ArgumentError
    If the number of design rows does not match the number of parameters.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function unit_design_to_physical(design::AbstractMatrix{<:Real}, params::Vector{ParameterBound})::Matrix{Float64}
    # Preliminary checks
    size(design, 1) == length(params) || throw(ArgumentError("design must have $(length(params)) rows"))

    # Convert each parameter using its own bounds and scale
    phys = similar(design, Float64)
    for (j, p) in enumerate(params)
        phys[j, :] .= from_unit(p, view(design, j, :); round_int = true)
    end
    return phys
end


"""
    ScreeningRecord

Outcome of the screening of one candidate.

# Fields
- `id::Int`: Column index in the design.
- `θ::Vector{Float64}`: Physical parameter values of the candidate.
- `status::Symbol`: `:accepted`, `:rejected` (error above threshold) or `:failed` (simulation failed).
- `error_V::Float64`: Worst error (V) over the datasets, according to `screening_norm` (`NaN` if failed).
- `dataset_errors_V::Vector{Float64}`: Error (V) on each screening dataset, in the order of `cfg.screening_datasets`.
"""
struct ScreeningRecord
    id::Int
    θ::Vector{Float64}
    status::Symbol
    error_V::Float64
    dataset_errors_V::Vector{Float64}
end


"""
    screening_error(U_sim, U_exp, norm::Symbol) -> Float64

Compute the voltage discrepancy between simulated and measured curves.

Parameters
----------
U_sim : AbstractVector{<:Real}
    Simulated cell voltages in volts.
U_exp : AbstractVector{<:Real}
    Measured cell voltages in volts.
norm : Symbol
    Error norm: `:linf` for maximum absolute deviation or `:rmse` for root-mean-square deviation.

Returns
-------
Float64
    Voltage error in volts, or `NaN` for empty or non-finite simulated values.

Raises
------
DimensionMismatch
    If the vectors have different lengths.
ArgumentError
    If `norm` is neither `:linf` nor `:rmse`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
Affonso Nobrega, marapendi, for the RMSE alternative.
"""
function screening_error(U_sim::AbstractVector{<:Real}, U_exp::AbstractVector{<:Real}, norm::Symbol)::Float64
    # Preliminary checks
    length(U_sim) == length(U_exp) || throw(DimensionMismatch("U_sim and U_exp must have the same length"))
    norm in (:linf, :rmse) || throw(ArgumentError("norm must be :linf or :rmse (got $norm)"))
    isempty(U_sim) && return NaN
    all(isfinite, U_sim) || return NaN

    # Compute the selected voltage-error norm
    d = U_sim .- U_exp
    if norm == :linf
        return Float64(maximum(abs, d))
    end
    return Float64(sqrt(sum(abs2, d) / length(d)))
end


"""
    screen_candidates(phys_design, datasets, cfg; simulate_fn) -> Vector{ScreeningRecord}

Screen every candidate against the experimental datasets. There is no stopping criterion: all
candidates are simulated and the number accepted is determined by the screening.

`simulate_fn(θ, dataset)` returns `(U_sim, U_exp)` (V, same length) or `nothing` / throws on failure.
A candidate is `:accepted` if its error is `≤ cfg.screening_threshold_V` on every dataset. Results
are stored per candidate index, so they do not depend on the number of threads.

Parameters
----------
phys_design : AbstractMatrix{<:Real}
    Physical parameter values, one candidate per column.
datasets : Vector{ScreeningDataset}
    Experimental datasets against which every candidate is evaluated.
cfg : ExtendedLocalSensitivityConfig
    Screening threshold, error norm, and execution settings.
simulate_fn : Function
    Callable `(θ, dataset) -> (U_sim, U_exp)`; may return `nothing` or throw when a simulation fails.

Returns
-------
Vector{ScreeningRecord}
    Screening status and voltage errors for each candidate, in design-column order.

Raises
------
ArgumentError
    If `datasets` is empty or `phys_design` contains no candidates.
ErrorException
    If fewer than 10% of candidates are accepted.
"""
function screen_candidates(phys_design::AbstractMatrix{<:Real},
                           datasets::Vector{ScreeningDataset},
                           cfg::ExtendedLocalSensitivityConfig;
                           simulate_fn::Function)::Vector{ScreeningRecord}
    # Preliminary checks
    !isempty(datasets) || throw(ArgumentError("datasets must not be empty"))
    n = size(phys_design, 2)
    n > 0 || throw(ArgumentError("phys_design must contain at least one candidate"))

    # Prepare candidate-indexed storage to preserve the design order when computed in parallel
    records = Vector{ScreeningRecord}(undef, n)

    function screen_one(k::Int)
        # Initialize this candidate's measurements and optimistic status
        θ = Vector{Float64}(phys_design[:, k]) # θ is this candidate's physical parameter vector.
        errs = fill(NaN, length(datasets))
        status = :accepted

        # Evaluate the candidate against each experimental dataset
        for (m, ds) in enumerate(datasets)
            res = try
                simulate_fn(θ, ds)
            catch
                nothing
            end
            e = res === nothing ? NaN : screening_error(res[1], res[2], cfg.screening_norm)
            errs[m] = e
            # Distinguish a failed simulation from a valid candidate above the threshold.
            if !isfinite(e)
                status = :failed
                break
            end
            e > cfg.screening_threshold_V && (status = :rejected)
        end
        worst = status == :failed ? NaN : maximum(errs)
        return ScreeningRecord(k, θ, status, worst, errs)
    end

    # Run serially or in parallel
    if cfg.parallel && Threads.nthreads() > 1
        # Wait for all candidates before computing the acceptance rate.
        @sync for k in 1:n
            Threads.@spawn records[k] = screen_one(k)
        end
    else
        for k in 1:n
            records[k] = screen_one(k)
        end
    end

    # Summarize statuses and check the acceptance rate
    n_acc = count(r -> r.status == :accepted, records)
    n_fail = count(r -> r.status == :failed, records)
    @info "Screening: $n_acc accepted, $(n - n_acc - n_fail) rejected, $n_fail failed (of $n candidates)."
    acceptance_rate = n_acc / n
    acceptance_rate < 0.1 &&
        error("Screening acceptance rate is $(round(100 * acceptance_rate; digits = 1)) %, below the 10 % minimum. Consider loosening screening_threshold_V or narrowing the parameter bounds (e.g. with ValidParameterRegion).")
    return records
end


"""
    alphapem_screening_simulator(params, cfg) -> Function

Construct a screening function that runs the real AlphaPEM model.

For each dataset, the stack's nominal physical parameters are overridden by `θ` for analysed
parameters only; all other values remain nominal. The polarisation curve is simulated at the
dataset's calibration currents.

Parameters
----------
params : Vector{ParameterBound}
    Analysed parameters, in the same order as the entries of `θ`.
cfg : ExtendedLocalSensitivityConfig
    Simulation settings, including gas-channel resolution and timeout.

Returns
-------
Function
    Callable `(θ, dataset) -> (U_sim, U_exp)` with voltages in volts, or `nothing` when simulation
    results cannot be extracted.
"""
function alphapem_screening_simulator(params::Vector{ParameterBound}, cfg::ExtendedLocalSensitivityConfig)
    names = [p.name for p in params]
    return function (θ::AbstractVector{<:Real}, ds::ScreeningDataset)
        # Load the dataset's nominal stack and experimental data
        fc = create_fuelcell(ds.fuel_cell_type, ds.voltage_zone; year = ds.year, nb_gc = cfg.nb_gc)
        exp_data = fc.pola_exp_data_cali

        # Override analysed parameters while keeping all others nominal
        overrides = (; (n => (p.type == :int ? Int(round(v)) : Float64(v))
                        for (n, p, v) in zip(names, params, θ))...)
        base = fc.physical_parameters
        base_nt = (; (f => getfield(base, f) for f in fieldnames(PhysicalParams))...)
        phys = PhysicalParams(; merge(base_nt, overrides)...)

        # Configure and run the simulation at the calibration currents
        sim_cfg = SimulationConfig(
            type_fuel_cell = ds.fuel_cell_type,
            year = ds.year,
            voltage_zone = ds.voltage_zone,
            type_current = PolarizationCalibrationParams(),
            numerical_parameters = NumericalParams(nb_gc = cfg.nb_gc,
                                                   max_run_time_s = cfg.max_run_time_s),
            type_display = :no_display,
            physical_parameters = phys,
        )
        simu = run_simulation(sim_cfg)
        simu.outputs === nothing && return nothing
        # Extract simulated voltages and align them with the measured points
        _, U_sim = _polarization_points_cali(simu.outputs, simu.current_density, exp_data.i_exp)
        length(U_sim) == length(exp_data.U_exp) || return nothing
        return (Float64.(U_sim), Float64.(exp_data.U_exp))
    end
end
