using AlphaPEM
using Test

const ELSA = AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

@testset "ExtendedLocalSensitivityAnalysis parameters" begin
    cfg = ELSA.ExtendedLocalSensitivityConfig()
    params, fixed = ELSA.select_parameters(cfg)

    # All PARAMETER_BOUNDS except the integer parameter :e
    @test length(params) == length(AlphaPEM.Config.PARAMETER_BOUNDS) - 1
    @test fixed == [:e]
    @test all(p -> p.min < p.max, params)
    @test all(p -> p.latex != "", params)
    @test any(p -> p.scale == :log, params)

    # Round-trip bijection, linear and log
    for p in params
        θ = [p.min, (p.min + p.max) / 2, p.max]
        @test ELSA.from_unit(p, ELSA.to_unit(p, θ)) ≈ θ
        @test ELSA.to_unit(p, θ[1]) ≈ 0.0 atol = 1e-12
        @test ELSA.to_unit(p, θ[3]) ≈ 1.0 atol = 1e-12
    end

    # :round policy keeps :e and rounds its levels
    cfg_round = ELSA.ExtendedLocalSensitivityConfig(integer_policy = :round)
    params_r, fixed_r = ELSA.select_parameters(cfg_round)
    @test isempty(fixed_r)
    pe = only(filter(p -> p.name == :e, params_r))
    @test ELSA.from_unit(pe, 0.4; round_int = true) == 4.0

    # Subset selection and invalid-input handling
    sub, _ = ELSA.select_parameters(ELSA.ExtendedLocalSensitivityConfig(parameter_subset = [:Hmem, :Re]))
    @test Set(p.name for p in sub) == Set([:Hmem, :Re])
    @test_throws ArgumentError ELSA.select_parameters(
        ELSA.ExtendedLocalSensitivityConfig(parameter_subset = [:nope]))
    bad = AlphaPEM.Parametrisation.Calibration.ParameterBound(:x, 0.0, 10.0, :real, :log, "", "")
    @test_throws ArgumentError ELSA.to_unit(bad, 1.0)
end

@testset "ExtendedLocalSensitivityAnalysis sensitivity (eq. 7)" begin
    θ_k = [0.3, 0.6, 0.9]

    @test ELSA.equidistant_levels(7) ≈ [0, 1/6, 2/6, 3/6, 4/6, 5/6, 1]

    # Linear model y = 2 + 3θ1 - θ2: exact slope normalised by the mean output over the sweep
    lin(θ) = [2 + 3θ[1] - θ[2]]
    S = ELSA.extended_local_sensitivity_matrix(lin, θ_k, 7)
    @test size(S) == (3, 1)
    @test S[1, 1] ≈ 3 / (2 + 3 * 0.5 - 0.6)
    @test S[2, 1] ≈ -1 / (2 + 3 * 0.3 - 0.5)
    @test S[3, 1] ≈ 0.0 atol = 1e-12

    # Several measurements
    multi(θ) = [1 + θ[1], 4 + 2θ[2]]
    S2 = ELSA.extended_local_sensitivity_matrix(multi, θ_k, 5)
    @test size(S2) == (3, 2)
    @test S2[1, 1] ≈ 1 / 1.5
    @test S2[2, 2] ≈ 2 / 5

    # Failed levels are removed; fewer than 2 valid levels gives NaN
    failing(θ) = θ[1] > 0.5 ? nothing : [2 + 3θ[1]]
    S3 = ELSA.extended_local_sensitivity_matrix(failing, [0.2], 7)
    @test S3[1, 1] ≈ 3 / (2 + 3 * 0.25)   # valid levels 0, 1/6, 2/6, 3/6; mean level 0.25
    onlyone(θ) = θ[1] < 0.01 ? [1.0] : nothing
    @test isnan(ELSA.extended_local_sensitivity_matrix(onlyone, [0.2], 7)[1, 1])

    @test_throws ArgumentError ELSA.extended_local_sensitivity_matrix(θ -> nothing, θ_k, 7)
    @test_throws ArgumentError ELSA.validate_config(ELSA.ExtendedLocalSensitivityConfig(n_levels = 2))
end

@testset "ExtendedLocalSensitivityAnalysis statistics (eq. 8-9)" begin
    # 3 points × 2 parameters × 2 measurements
    S_n = zeros(3, 2, 2)
    S_n[:, 1, :] = [3.0 4.0; 6.0 8.0; 0.0 10.0]   # norms 5, 10, 10
    S_n[:, 2, :] = [1.0 0.0; 1.0 0.0; 1.0 0.0]    # norms 1, 1, 1
    M = ELSA.output_sensitivity_magnitude(S_n)
    @test size(M) == (3, 2)
    @test M[:, 1] ≈ [5.0, 10.0, 10.0]
    @test M[:, 2] ≈ [1.0, 1.0, 1.0]

    s = ELSA.sensitivity_statistics(S_n)
    @test s.median ≈ [10.0, 1.0]
    @test s.std[2] ≈ 0.0 atol = 1e-12
    @test s.std[1] ≈ sqrt(sum((M[:, 1] .- 25 / 3) .^ 2) / 2)

    # Failed points (NaN) are ignored
    S_nan = copy(S_n)
    S_nan[1, 1, :] .= NaN
    s2 = ELSA.sensitivity_statistics(S_nan)
    @test s2.median ≈ [10.0, 1.0]
    @test s2.std[1] ≈ 0.0 atol = 1e-12   # two equal magnitudes (10, 10)
    S_nan[:, 2, :] .= NaN
    @test isnan(ELSA.sensitivity_statistics(S_nan).median[2])
end

@testset "ExtendedLocalSensitivityAnalysis collinearity (eq. 10)" begin
    # Rows: 1 and 2 are collinear, 3 is orthogonal, and 4 is null
    S = [1.0 2.0 0.0;
         3.0 6.0 0.0;
         0.0 0.0 5.0;
         0.0 0.0 0.0]
    Ψ = ELSA.collinearity_index(S)
    @test size(Ψ) == (4, 4)
    @test Ψ[1, 2] ≈ 1.0
    @test Ψ[1, 3] ≈ 0.0 atol = 1e-12
    @test Ψ ≈ Ψ'
    @test Ψ[1, 1] ≈ 1.0
    @test all(Ψ[4, :] .== 0.0)                      # null parameter: floor avoids 0/0
    @test ELSA.collinearity_index([1.0 -2.0; -2.0 4.0])[1, 2] ≈ 1.0
    @test 0.0 <= ELSA.collinearity_index([1.0 2.0 3.0; 2.0 -1.0 0.5])[1, 2] <= 1.0

    # Per-point indices and their median
    S_n = zeros(3, 2, 2)
    S_n[1, :, :] = [1.0 0.0; 1.0 0.0]               # ψ = 1
    S_n[2, :, :] = [1.0 0.0; 0.0 1.0]               # ψ = 0
    S_n[3, :, :] = [1.0 0.0; 1.0 0.0]               # ψ = 1
    Ψs = ELSA.collinearity_indices(S_n)
    @test size(Ψs) == (3, 2, 2)
    @test Ψs[:, 1, 2] ≈ [1.0, 0.0, 1.0]
    @test ELSA.median_collinearity(S_n)[1, 2] ≈ 1.0

    # Failed points are ignored when taking the median
    S_nan = copy(S_n)
    S_nan[1, 1, :] .= NaN
    @test isnan(ELSA.collinearity_indices(S_nan)[1, 1, 2])
    @test ELSA.median_collinearity(S_nan)[1, 2] ≈ 0.5

    @test ELSA.influential_mask([0.5, 0.05, NaN, 0.1], 0.1) == BitVector([true, false, false, true])
    @test ELSA.influential_mask([Inf, -Inf, 0.2], 0.1) == BitVector([false, false, true])
    @test_throws ArgumentError ELSA.influential_mask([0.2], -0.1)
    @test_throws ArgumentError ELSA.influential_mask([0.2], Inf)
end

@testset "ExtendedLocalSensitivityAnalysis sampling and screening (step 4)" begin
    @testset "scrambled Sobol design" begin
        D1 = ELSA.generate_unit_design(64, 5; seed = 1)
        D2 = ELSA.generate_unit_design(64, 5; seed = 1)
        D3 = ELSA.generate_unit_design(64, 5; seed = 2)
        @test size(D1) == (5, 64)
        @test D1 == D2
        @test D1 != D3
        @test all(0 .<= D1 .<= 1)
        @test !any(all(D1[:, k] .== 0) for k in 1:64)   # never the all-lower-bound point
        @test size(ELSA.generate_unit_design(2^9, 4; seed = 1)) == (4, 2^9)
        @test size(ELSA.generate_unit_design(350, 4; seed = 1)) == (4, 2^9)
        D4 = ELSA.generate_unit_design(8, 3; seed = 1, nominal_unit = [0.5, 0.5, 0.5])
        @test size(D4) == (3, 9)
        @test D4[:, 1] == [0.5, 0.5, 0.5]
        @test_throws ArgumentError ELSA.generate_unit_design(8, 3; seed = 1, nominal_unit = [0.5])
    end

    @testset "physical mapping" begin
        params = ELSA.select_parameters(ELSA.ExtendedLocalSensitivityConfig())[1][1:3]
        D = ELSA.generate_unit_design(16, 3; seed = 3)
        P = ELSA.unit_design_to_physical(D, params)
        for (j, p) in enumerate(params)
            @test all(p.min - 1e-9 .<= P[j, :] .<= p.max + 1e-9)
        end
    end

    @testset "screening" begin
        @test ELSA.screening_error([1.0, 2.0], [1.1, 2.3], :linf) ≈ 0.3
        @test ELSA.screening_error([1.0, 2.0], [1.0, 2.2], :rmse) ≈ sqrt(0.02)
        @test isnan(ELSA.screening_error([NaN, 1.0], [1.0, 1.0], :linf))

        ds = [ELSA.ScreeningDataset(:ZSW_nominal, 2024, :full),
              ELSA.ScreeningDataset(:ZSW_Pa_1_61_Pc_1_41, 2024, :full)]
        # Candidate k has a constant bias θ[1] on every dataset; negative values simulate failure.
        phys = [0.1 0.5 -1.0 0.2; 0.0 0.0 0.0 0.0]
        sim(θ, d) = θ[1] < 0 ? nothing : (fill(0.7 + θ[1], 3), fill(0.7, 3))
        for parallel in (false, true)
            cfg = ELSA.ExtendedLocalSensitivityConfig(parallel = parallel, screening_datasets = ds,
                                                      screening_threshold_V = 0.3)
            rec = ELSA.screen_candidates(phys, ds, cfg; simulate_fn = sim)
            @test [r.status for r in rec] == [:accepted, :rejected, :failed, :accepted]
            @test rec[1].error_V ≈ 0.1
            @test length(rec[1].dataset_errors_V) == 2
            @test [r.id for r in rec] == 1:4
        end
        cfg = ELSA.ExtendedLocalSensitivityConfig(screening_datasets = ds)
        @test_throws ErrorException ELSA.screen_candidates(phys[:, 2:3], ds, cfg; simulate_fn = sim)

        # Exactly 10% acceptance is allowed; a lower rate aborts screening.
        threshold_design = zeros(1, 10)
        pass_first(θ, d) = θ[1] == 1.0 ? (fill(0.7, 3), fill(0.7, 3)) :
                                         (fill(1.1, 3), fill(0.7, 3))
        threshold_design[1, 1] = 1.0
        @test length(ELSA.screen_candidates(threshold_design, ds, cfg; simulate_fn = pass_first)) == 10
        below_threshold_design = hcat(threshold_design, zeros(1, 1))
        @test_throws ErrorException ELSA.screen_candidates(
            below_threshold_design, ds, cfg; simulate_fn = pass_first)
        @test_throws ArgumentError ELSA.screen_candidates(zeros(1, 0), ds, cfg; simulate_fn = sim)
        @test_throws ArgumentError ELSA.screen_candidates(zeros(1, 1), ELSA.ScreeningDataset[], cfg;
                                                          simulate_fn = sim)
    end

    @testset "configuration" begin
        cfg = ELSA.ExtendedLocalSensitivityConfig()
        @test cfg.n_candidates == 2^9
        @test cfg.include_nominal == false
        @test !hasproperty(cfg, :scramble) && !hasproperty(cfg, :n_target)
        @test ELSA.validate_config(ELSA.ExtendedLocalSensitivityConfig(n_candidates = 350)).n_candidates == 350
        @test_throws Exception ELSA.validate_config(
            ELSA.ExtendedLocalSensitivityConfig(screening_datasets = ELSA.ScreeningDataset[]))
    end
end

@testset "ExtendedLocalSensitivityAnalysis operating cases (step 5)" begin
    default_bounds = copy(AlphaPEM.Config.OPERATING_CONDITIONS_BOUNDS)
    cfg = ELSA.ExtendedLocalSensitivityConfig()
    cases = ELSA.generate_operating_case_library(cfg)

    @test length(cases) == 3^length(ELSA.OPERATING_CONDITIONS_NAMES_FOR_PARAMETRISATION) == 729
    @test first(cases).label == "factorial_001"
    @test last(cases).label == "factorial_729"
    @test length(Set(case.label for case in cases)) == 729
    @test all(case -> case.operating_conditions.Pa_des -
                      case.operating_conditions.Pc_des ≈ 0.2e5, cases)
    nominal_base = AlphaPEM.Fuelcell.create_fuelcell(
        cfg.fuel_cell_type, cfg.voltage_zone; year = cfg.year, nb_gc = cfg.nb_gc).operating_conditions
    @test all(case -> case.operating_conditions.i_min_stoich == nominal_base.i_min_stoich, cases)
    @test AlphaPEM.Config.OPERATING_CONDITIONS_BOUNDS == default_bounds
    @test cfg.current_density_grid == ELSA.DEFAULT_CURRENT_DENSITY_GRID

    # Each parameterized operating condition spans its effective minimum, midpoint, and maximum.
    for name in ELSA.OPERATING_CONDITIONS_NAMES_FOR_PARAMETRISATION
        lo, hi, _ = AlphaPEM.Config.OPERATING_CONDITIONS_BOUNDS[name]
        values = Set(getfield(case.operating_conditions, name) for case in cases)
        @test values == Set([lo, (lo + hi) / 2.0, hi])
    end

    # Explicit base values are retained for operating conditions not included in the parameterization.
    base = AlphaPEM.Config.OperatingConditions(i_min_stoich = 0.7, y_H2_in = 0.9)
    cfg_custom = ELSA.ExtendedLocalSensitivityConfig(
        base_operating_conditions = base,
        pressure_gap = 0.1e5,
        operating_bounds_overrides = Dict(:T_des => (300.0, 340.0, :real)),
    )
    custom_cases = ELSA.generate_operating_case_library(cfg_custom)
    @test all(case -> case.operating_conditions.i_min_stoich == 0.7, custom_cases)
    @test all(case -> case.operating_conditions.y_H2_in == 0.9, custom_cases)
    @test all(case -> case.operating_conditions.Pa_des -
                      case.operating_conditions.Pc_des ≈ 0.1e5, custom_cases)
    @test Set(case.operating_conditions.T_des for case in custom_cases) ==
          Set([300.0, 320.0, 340.0])

    # Invalid pressure-derived values are rejected rather than clamped.
    cfg_bad_pc_bounds = ELSA.ExtendedLocalSensitivityConfig(
        operating_bounds_overrides = Dict(:Pc_des => (1.4e5, 3.5e5, :real)),
    )
    @test_throws ArgumentError ELSA.generate_operating_case_library(cfg_bad_pc_bounds)
    @test_throws ArgumentError ELSA.validate_config(
        ELSA.ExtendedLocalSensitivityConfig(current_density_grid = [0.0, 0.2, 0.1]))
    @test_throws ArgumentError ELSA.validate_config(
        ELSA.ExtendedLocalSensitivityConfig(operating_bounds_overrides = Dict(:unknown => (0.0, 1.0, :real))))

    pressure_constraint = AlphaPEM.Config.exact_pressure_gap_constraint(0.2e5)
    @test pressure_constraint.target == :Pc_des
    @test pressure_constraint.kind == :(=)
    @test pressure_constraint.fn(Dict(:Pa_des => 2.0e5)) ≈ 1.8e5
    @test length(AlphaPEM.Config.default_operating_condition_constraints()) == 2
    @test_throws ArgumentError AlphaPEM.Config.exact_pressure_gap_constraint(0.0)
    @test_throws ArgumentError AlphaPEM.Config.exact_pressure_gap_constraint(Inf)
end

@testset "ExtendedLocalSensitivityAnalysis ascending current profile (step 6.1)" begin
    AscendingParams = AlphaPEM.Config.AscendingPolarizationParams
    AscendingCurrent = AlphaPEM.Currents.AscendingPolarizationCurrent
    current = AlphaPEM.Currents.current

    # Initialize at 1 A/cm², return to OCV, then make one ascending pass.
    p = AscendingParams(delta_t_ini = 20.0, v_load = 1e4, delta_t_break = 10.0,
                        delta_t_break_lower_0_3 = 8.0, delta_t_break_OCV = 6.0,
                        delta_t_measurement = 5.0,
                        i_grid = [0.0, 0.2e4, 1.0e4, 2.0e4])
    c = AscendingCurrent(p)
    @test c.i_targets == p.i_grid
    @test length(c.di_transitions) == 5
    @test c.di_transitions[1:2] == [1.0e4, -1.0e4] # initial rise to 1 A/cm², then return to OCV
    @test all(c.di_transitions[3:end] .> 0.0)       # requested sweep only ascends
    @test sum(c.di_transitions) ≈ 2.0e4
    @test c.t_meas_ends[1] ≈ 1.0 + 20.0 + 1.0 + 6.0 # initial hold, return ramp and OCV hold
    @test c.t_meas_ends[2] - c.t_meas_ends[1] ≈ 0.2 + 8.0    # ramp + stabilisation up to 0.3 A/cm²
    @test c.t_meas_ends[3] - c.t_meas_ends[2] ≈ 0.8 + 10.0   # ramp + stabilisation above 0.3 A/cm²
    @test c.time_interval == (0.0, c.t_meas_ends[end])
    @test all(c.t_meas_ends .- c.t_meas_starts .≈ 5.0)
    for k in eachindex(c.i_targets)
        @test current(c, c.t_meas_starts[k]) ≈ c.i_targets[k] atol = 1e-6
        # The next ramp begins at the window end, with a leading tanh tail below 0.1 % of its step.
        @test current(c, c.t_meas_ends[k]) ≈ c.i_targets[k] atol = 10.0
    end

    # Even when OCV is not a requested data point, initialize at 1 A/cm² and return to OCV first.
    c_ramp = AscendingCurrent(AscendingParams(
        delta_t_ini = 20.0, v_load = 1e4, delta_t_break = 10.0, delta_t_break_lower_0_3 = 10.0,
        delta_t_break_OCV = 10.0, delta_t_measurement = 5.0, i_grid = [0.05e4, 1.0e4]))
    @test length(c_ramp.di_transitions) == 4
    @test c_ramp.t_starts[1:2] == [0.0, 21.0]
    @test c_ramp.t_meas_ends[1] ≈ 1.0 + 20.0 + 1.0 + 10.0 + 0.05 + 10.0

    # Invalid profiles are rejected
    @test_throws ArgumentError AscendingCurrent(AscendingParams())
    @test_throws ArgumentError AscendingCurrent(AscendingParams(i_grid = [1e4, 5e3]))
    @test_throws ArgumentError AscendingCurrent(AscendingParams(i_grid = [-1.0, 5e3]))
    @test_throws ArgumentError AscendingCurrent(AscendingParams(i_grid = [0.0, 1e4], v_load = 0.0))
    @test_throws ArgumentError AscendingCurrent(AscendingParams(
        i_grid = [0.0, 1e4], delta_t_break = 2.0, delta_t_measurement = 5.0))
    @test_throws ArgumentError AscendingCurrent(AscendingParams(
        i_grid = [0.0, 1e4], delta_t_break_OCV = 2.0, delta_t_measurement = 5.0))
    @test AlphaPEM.Currents.create_current(p) isa AscendingCurrent
    sim_cfg = AlphaPEM.Config.SimulationConfig(type_current = p)
    @test AlphaPEM.Config.validate_config(sim_cfg) === sim_cfg

    # The run configuration keeps AlphaPEM's durations and converts the grid from A/cm² to A/m²
    cfg = ELSA.ExtendedLocalSensitivityConfig()
    params = ELSA.ascending_polarization_params(cfg)
    cali = AlphaPEM.Config.PolarizationCalibrationParams()
    @test params.i_grid ≈ cfg.current_density_grid .* 1e4
    @test (params.delta_t_ini, params.v_load, params.delta_t_break,
           params.delta_t_break_lower_0_3, params.delta_t_break_OCV,
           params.delta_t_measurement) ==
          (cali.delta_t_ini, cali.v_load, cali.delta_t_break,
           cali.delta_t_break_lower_0_3, cali.delta_t_break_OCV,
           cali.delta_t_measurement)
    @test length(AscendingCurrent(params).i_targets) == 15
    @test_throws ArgumentError ELSA.validate_config(ELSA.ExtendedLocalSensitivityConfig(v_load = 0.0))
    @test_throws ArgumentError ELSA.validate_config(
        ELSA.ExtendedLocalSensitivityConfig(delta_t_measurement = 400.0))
    @test_throws ArgumentError ELSA.validate_config(
        ELSA.ExtendedLocalSensitivityConfig(delta_t_break_OCV = -1.0))
    @test_throws ArgumentError ELSA.validate_config(
        ELSA.ExtendedLocalSensitivityConfig(delta_t_ini = NaN))
end

@testset "ExtendedLocalSensitivityAnalysis voltage curve (step 6.1)" begin
    AscendingParams = AlphaPEM.Config.AscendingPolarizationParams
    c = AlphaPEM.Currents.AscendingPolarizationCurrent(AscendingParams(
        delta_t_ini = 20.0, v_load = 1e4, delta_t_break = 10.0, delta_t_break_lower_0_3 = 8.0,
        delta_t_break_OCV = 6.0,
        delta_t_measurement = 5.0, i_grid = [0.0, 0.2e4, 1.0e4, 2.0e4]))
    t = collect(0.0:0.05:c.time_interval[2])
    U = 0.95 .- 1e-5 .* AlphaPEM.Currents.current.(Ref(c), t)
    expected_U = 0.95 .- 1e-5 .* c.i_targets

    # Complete curve: voltages at the requested current densities, without interpolation
    full = ELSA.voltage_curve_from_history(t, U, c)
    @test full.status == :complete
    @test all(full.valid)
    @test full.current_density == c.i_targets ./ 1e4
    @test full.voltage ≈ expected_U atol = 1e-5   # the last window sample sees the next ramp's tanh tail
    @test full.last_reached_current_density == 2.0

    # Solver stopped during the third stabilisation: the first two points are kept
    stopped = t .<= c.t_meas_starts[3] + 1.0
    partial = ELSA.voltage_curve_from_history(t[stopped], U[stopped], c)
    @test partial.status == :partial
    @test partial.valid == [true, true, false, false]
    @test partial.voltage[1:2] ≈ expected_U[1:2] atol = 1e-5
    @test all(isnan, partial.voltage[3:4])
    @test partial.last_reached_current_density == 0.2

    # Solver stopped before the first measurement window ends, or no output at all
    early = t .<= c.t_meas_ends[1] - 1.0
    @test ELSA.voltage_curve_from_history(t[early], U[early], c).status == :invalid
    no_points = ELSA.voltage_curve_from_history(Float64[], Float64[], c)
    @test no_points.status == :invalid
    @test isnan(no_points.last_reached_current_density)

    # A non-finite voltage ends the valid prefix; later points are never recovered
    U_bad_2 = copy(U)
    U_bad_2[findfirst(t .>= c.t_meas_starts[2] + 1.0)] = NaN
    @test ELSA.voltage_curve_from_history(t, U_bad_2, c).valid == [true, false, false, false]
    U_bad_3 = copy(U)
    U_bad_3[findfirst(t .>= c.t_meas_starts[3] + 1.0)] = Inf
    @test ELSA.voltage_curve_from_history(t, U_bad_3, c).valid == [true, true, false, false]
    U_bad_1 = copy(U)
    U_bad_1[findfirst(t .>= c.t_meas_starts[1] + 1.0)] = NaN
    @test ELSA.voltage_curve_from_history(t, U_bad_1, c).status == :invalid

    # Invalid inputs
    @test_throws ArgumentError ELSA.voltage_curve_from_history(t, U[1:end-1], c)
    @test_throws ArgumentError ELSA.VoltageCurve([1.0], [1.0, 2.0], trues(1))

    # Simulators without outputs give an invalid curve; other profiles are refused
    no_outputs = ELSA.extract_voltage_curve((current_density = c, outputs = nothing))
    @test no_outputs.status == :invalid
    @test length(no_outputs.voltage) == 4
    @test_throws ArgumentError ELSA.extract_voltage_curve((current_density = 1.0, outputs = nothing))
end

@testset "ExtendedLocalSensitivityAnalysis ascending sweep on AlphaPEM (step 6.1)" begin
    numerical_params = AlphaPEM.Config.NumericalParams(
        nb_gc = 5, nb_gdl = 2, nb_mpl = 1, nb_man = 1, rtol = 1e-2, atol = 1e-4,
        maxiters = 1e4, max_run_time_s = 60.0, save_freq = 100.0)
    run_sweep(grid; v_load) = begin
        cfg = ELSA.ExtendedLocalSensitivityConfig(
            current_density_grid = grid, delta_t_ini = 10.0, v_load = v_load, delta_t_break = 10.0,
            delta_t_break_lower_0_3 = 10.0, delta_t_break_OCV = 10.0, delta_t_measurement = 5.0)
        sim_cfg = AlphaPEM.Config.SimulationConfig(
            type_fuel_cell = :ZSW_nominal,
            type_current = ELSA.ascending_polarization_params(cfg),
            numerical_parameters = numerical_params,
            voltage_zone = :full,
            operating_conditions = AlphaPEM.Config.OperatingConditions(Sa = 3.0, Sc = 4.0),
            type_display = :no_display)
        return AlphaPEM.Application.run_simulation(sim_cfg)
    end

    # Reachable currents: complete curve with a decreasing polarization voltage
    complete = ELSA.extract_voltage_curve(run_sweep([0.0, 0.05, 0.1]; v_load = 0.01e4))
    @test complete.status == :complete
    @test issorted(complete.voltage; rev = true)
    @test all(0.5 .< complete.voltage .< 1.1)

    # Currents far beyond the cell capability: the solver fails, the reached prefix is kept
    stopped = ELSA.extract_voltage_curve(run_sweep([0.0, 0.5, 4.0, 8.0]; v_load = 0.1e4))
    @test stopped.status == :partial
    @test stopped.valid[1]
    @test !stopped.valid[end]
    @test issorted(stopped.valid; rev = true)
    @test all(isnan, stopped.voltage[.!stopped.valid])
end

@testset "ExtendedLocalSensitivityAnalysis current-density summaries (step 6.2)" begin
    bounds = AlphaPEM.Parametrisation.Calibration.ParameterBound
    params = [
        bounds(:Hmem, 0.1, 10.0, :real, :log, "Hmem", "Membrane thickness"),
        bounds(:Re, 0.1, 10.0, :real, :log, "Re", "Electrical resistance"),
    ]
    grid = [0.0, 0.25, 0.5, 0.75, 1.5, 1.6, 2.0]
    S = [
        -1.0  1.0  -2.0  3.0  NaN  -4.0  NaN
        NaN   NaN   NaN   NaN  NaN   NaN  NaN
    ]
    original_S = copy(S)
    cfg = ELSA.ExtendedLocalSensitivityConfig(current_density_grid = grid)

    summary = ELSA.summarize_current_density_sensitivities(S, params, cfg)
    @test size(summary, 1) == length(params) * (1 + length(cfg.current_density_zones))
    @test Set(summary.region) == Set([:all_current_densities; [z.name for z in cfg.current_density_zones]])
    @test isequal(S, original_S) # Keep each positive/negative voltage response unchanged in the raw matrix.

    row(parameter, region) = summary[(summary.parameter .== parameter) .&
                                     (summary.region .== region), :]
    global_hmem = row(:Hmem, :all_current_densities)
    @test only(global_hmem.median_magnitude) == 2.0
    @test only(global_hmem.median_signed) == -1.0
    @test only(global_hmem.valid_points) == 5
    @test only(global_hmem.total_points) == 7
    @test only(global_hmem.coverage) ≈ 5 / 7

    low = row(:Hmem, :low_current_density_or_kinetic)
    @test only(low.median_magnitude) == 1.0
    @test only(low.median_signed) == -1.0
    @test only(low.valid_points) == 3
    @test only(low.total_points) == 3
    @test only(low.coverage) == 1.0

    medium = row(:Hmem, :medium_current_density_or_ohmic)
    @test only(medium.median_magnitude) == 3.0
    @test only(medium.median_signed) == 3.0
    @test only(medium.valid_points) == 1
    @test only(medium.total_points) == 2
    @test only(medium.coverage) == 0.5

    high = row(:Hmem, :high_current_density_or_mass_transport)
    @test only(high.median_magnitude) == 4.0
    @test only(high.median_signed) == -4.0
    @test only(high.valid_points) == 1
    @test only(high.total_points) == 2
    @test only(high.coverage) == 0.5

    empty_region = row(:Re, :all_current_densities)
    @test isnan(only(empty_region.median_magnitude))
    @test isnan(only(empty_region.median_signed))
    @test only(empty_region.valid_points) == 0
    @test only(empty_region.coverage) == 0.0

    no_grid_points_cfg = ELSA.ExtendedLocalSensitivityConfig(
        current_density_grid = grid,
        current_density_zones = [ELSA.CurrentDensityZone(:above_grid, 3.0, Inf)])
    no_grid_points = ELSA.summarize_current_density_sensitivities(S, params, no_grid_points_cfg)
    above_grid = no_grid_points[(no_grid_points.parameter .== :Hmem) .&
                                (no_grid_points.region .== :above_grid), :]
    @test only(above_grid.valid_points) == 0
    @test only(above_grid.total_points) == 0
    @test isnan(only(above_grid.coverage))
    @test isnan(only(above_grid.median_magnitude))

    # Zones may be customized, but overlapping intervals and duplicate names are rejected.
    custom_zones = [ELSA.CurrentDensityZone(:first, 0.0, 1.0),
                    ELSA.CurrentDensityZone(:second, 0.5, 2.0)]
    @test_throws ArgumentError ELSA.validate_config(
        ELSA.ExtendedLocalSensitivityConfig(current_density_zones = custom_zones))
    duplicate_zones = [ELSA.CurrentDensityZone(:same, 0.0, 1.0),
                       ELSA.CurrentDensityZone(:same, 1.0, 2.0)]
    @test_throws ArgumentError ELSA.validate_config(
        ELSA.ExtendedLocalSensitivityConfig(current_density_zones = duplicate_zones))

    @test_throws ArgumentError ELSA.summarize_current_density_sensitivities(
        zeros(1, length(grid)), params, cfg)
    @test_throws ArgumentError ELSA.summarize_current_density_sensitivities(
        zeros(length(params), length(grid) - 1), params, cfg)
end
