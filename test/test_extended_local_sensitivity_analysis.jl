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
