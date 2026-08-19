using AlphaPEM
using Test

# Helper: minimal numerical parameters for fast smoke tests.
function quick_numerical_params()
    return AlphaPEM.Config.NumericalParams(
        nb_gc = 1,
        nb_gdl = 2,
        nb_mpl = 1,
        nb_man = 1,
        rtol = 1e-2,
        atol = 1e-4,
        maxiters = 1e4,
        max_run_time_s = 30.0,
        save_freq = 100.0,
    )
end

@testset "AlphaPEM.jl" begin
    # Basic loading test
    @test true

    # Verify access to core modules and types
    @test isdefined(AlphaPEM, :Config)
    @test isdefined(AlphaPEM.Config, :SimulationConfig)
    @test isdefined(AlphaPEM, :Core)
    @test isdefined(AlphaPEM, :Application)
    @test isdefined(AlphaPEM, :Fuelcell)
    @test isdefined(AlphaPEM, :Currents)
    @test isdefined(AlphaPEM, :Parametrisation)
    @test isdefined(AlphaPEM, :Utils)
end

@testset "Config" begin
    using AlphaPEM.Config

    # Parameter structs can be instantiated with defaults
    np = NumericalParams()
    @test np.nb_gc >= 1

    sp = StepParams()
    @test sp.i_step > 0.0

    pp = PolarizationParams()
    @test pp.di_step > 0.0

    eis = EISParams()
    @test eis.nb_f > 0

    oc = OperatingConditions()
    @test oc.Sa > 0.0
    @test oc.Sc > 0.0

    # SimulationConfig default construction and validation
    cfg = SimulationConfig()
    @test cfg.type_fuel_cell == :ZSW_nominal
    @test validate_config(cfg) === cfg

    # Custom configuration
    cfg_custom = SimulationConfig(
        type_fuel_cell = :ZSW_nominal,
        type_current = StepParams(),
        numerical_parameters = quick_numerical_params(),
        type_display = :no_display,
    )
    @test cfg_custom.type_display == :no_display
    @test validate_config(cfg_custom) === cfg_custom

    # Invalid configuration should error
    @test_throws ErrorException validate_config(SimulationConfig(type_purge = :unknown_purge))
end

@testset "Fuelcell" begin
    using AlphaPEM.Fuelcell

    # Factory for the default ZSW stack
    fc = create_fuelcell(:ZSW_nominal, :full)
    @test fc isa AbstractFuelCell
    @test fc.physical_parameters isa AlphaPEM.Config.PhysicalParams
    @test fc.operating_conditions isa AlphaPEM.Config.OperatingConditions

    # Other supported fuel cell types
    @test create_fuelcell(:EH_nominal, :full) isa AbstractFuelCell
    @test create_fuelcell(:ZSW_T_84, :before_voltage_drop) isa AbstractFuelCell

    # Accessor functions
    @test physical_parameters(fc) isa AlphaPEM.Config.PhysicalParams
    @test operating_conditions(fc, :ZSW_nominal) isa AlphaPEM.Config.OperatingConditions
    @test undetermined_parameters(fc, :full) isa Vector{Tuple{Symbol, Float64, Float64}}
    @test !isempty(undetermined_parameters(fc, :full))
end

@testset "Currents" begin
    using AlphaPEM.Currents
    using AlphaPEM.Fuelcell: create_fuelcell

    fc = create_fuelcell(:ZSW_nominal, :full)

    # Factory for each current profile type
    step_current = create_current(StepParams(), fc)
    @test step_current isa StepCurrent

    pola_current = create_current(PolarizationParams(i_max = 1.5e4), fc)
    @test pola_current isa PolarizationCurrent

    cali_current = create_current(PolarizationCalibrationParams(), fc)
    @test cali_current isa PolarizationCalibrationCurrent

    eis_current = create_current(EISParams(nb_f = 3), fc)
    @test eis_current isa EISCurrent

    # Default factory dispatch
    @test create_current(StepParams()) isa StepCurrent
end

@testset "Simulation - step current" begin
    using AlphaPEM.Config: SimulationConfig, StepParams
    using AlphaPEM.Application: run_simulation

    cfg = SimulationConfig(
        type_fuel_cell = :ZSW_nominal,
        type_current = StepParams(
            delta_t_ini = 10.0,
            delta_t_load = 1.0,
            delta_t_break = 5.0,
            i_ini = 0.5e4,
            i_step = 1.0e4,
        ),
        numerical_parameters = quick_numerical_params(),
        type_display = :no_display,
        display_timing = :postrun,
    )

    simu = run_simulation(cfg)
    @test simu isa AlphaPEM.Core.Models.AlphaPEM
    @test simu.outputs !== nothing
    @test !isempty(simu.outputs.solver.t)
    @test !isempty(simu.outputs.solver.states)
end

@testset "Simulation - polarization curve" begin
    using AlphaPEM.Config: SimulationConfig, PolarizationParams
    using AlphaPEM.Application: run_simulation

    cfg = SimulationConfig(
        type_fuel_cell = :ZSW_nominal,
        type_current = PolarizationParams(
            delta_t_ini = 10.0,
            di_step = 0.5e4,
            v_load = 0.01e4,
            delta_t_break = 10.0,
            delta_t_break_lower_0_3 = 10.0,
            delta_t_break_OCV = 10.0,
            delta_t_measurement = 10.0,
            i_max = 1.5e4,
        ),
        numerical_parameters = quick_numerical_params(),
        voltage_zone = :full,
        type_display = :no_display,
        display_timing = :postrun,
    )

    simu = run_simulation(cfg)
    @test simu isa AlphaPEM.Core.Models.AlphaPEM
    @test simu.outputs !== nothing
    @test !isempty(simu.outputs.solver.t)
end

@testset "Simulation - EIS" begin
    using AlphaPEM.Config: SimulationConfig, EISParams
    using AlphaPEM.Application: run_simulation

    cfg = SimulationConfig(
        type_fuel_cell = :ZSW_nominal,
        type_current = EISParams(
            i_EIS = 0.5e4,
            ratio = 5.0 / 100.0,
            f_power_min = 0.0,
            f_power_max = 2.0,
            nb_f = 3,
            nb_points = 5,
        ),
        numerical_parameters = quick_numerical_params(),
        type_display = :no_display,
        display_timing = :postrun,
    )

    simu = run_simulation(cfg)
    @test simu isa AlphaPEM.Core.Models.AlphaPEM
    @test simu.outputs !== nothing
    @test !isempty(simu.outputs.solver.t)
end

@testset "Simulation - polarization calibration" begin
    using AlphaPEM.Config: SimulationConfig, PolarizationCalibrationParams
    using AlphaPEM.Application: run_simulation

    cfg = SimulationConfig(
        type_fuel_cell = :ZSW_nominal,
        type_current = PolarizationCalibrationParams(
            delta_t_ini = 10.0,
            v_load = 0.01e4,
            delta_t_break = 10.0,
            delta_t_break_lower_0_3 = 5.0,
            delta_t_break_OCV = 5.0,
            delta_t_measurement = 5.0,
        ),
        numerical_parameters = quick_numerical_params(),
        voltage_zone = :full,
        type_display = :no_display,
        display_timing = :postrun,
    )

    simu = run_simulation(cfg)
    @test simu isa AlphaPEM.Core.Models.AlphaPEM
    @test simu.outputs !== nothing
    @test !isempty(simu.outputs.solver.t)
end

@testset "Simulation - batch run" begin
    using AlphaPEM.Config: SimulationConfig, PolarizationParams
    using AlphaPEM.Application: run_simulation

    pola_params = PolarizationParams(
        delta_t_ini = 10.0,
        di_step = 0.5e4,
        v_load = 0.01e4,
        delta_t_break = 10.0,
        delta_t_break_lower_0_3 = 10.0,
        delta_t_break_OCV = 10.0,
        delta_t_measurement = 10.0,
        i_max = 1.5e4,
    )
    num_params = quick_numerical_params()

    cfgs = [
        SimulationConfig(
            type_fuel_cell = :ZSW_nominal,
            type_current = pola_params,
            numerical_parameters = num_params,
            voltage_zone = :full,
            type_display = :no_display,
            display_timing = :postrun,
        ),
        SimulationConfig(
            type_fuel_cell = :ZSW_T_84,
            type_current = pola_params,
            numerical_parameters = num_params,
            voltage_zone = :before_voltage_drop,
            type_display = :no_display,
            display_timing = :postrun,
        ),
    ]

    simus = run_simulation(cfgs)
    @test length(simus) == 2
    @test all(s -> s isa AlphaPEM.Core.Models.AlphaPEM, simus)
    @test all(s -> s.outputs !== nothing, simus)
end

@testset "Parametrisation.ValidParameterRegion" begin
    using AlphaPEM.Parametrisation.ValidParameterRegion
    using AlphaPEM.Config: ValidityCriteriaConfig

    # Sampling and bounds
    cfg = ValidityAnalysisConfig(
        fuel_cell_type = :ZSW_nominal,
        year = 2024,
        voltage_zone = :full,
        n_samples = 4,
        parallel = false,
        save_curves = false,
        hyperbox_finder_method = nothing,
        max_run_time_s = 30.0,
    )

    X, pb = generate_test_samples(cfg)
    @test size(X, 1) == cfg.n_samples
    @test size(X, 2) == length(pb.bounds)

    # Validity criteria on synthetic curves
    criteria = AlphaPEM.Parametrisation.ValidParameterRegion.ValidityCriteriaConfig()

    # Example from the module docstring: monotonically decreasing polarization curve
    valid_result = classify_polarization_curve(
        [1.05, 0.92, 0.80, 0.65, 0.50],
        [0.0, 5e3, 1e4, 1.5e4, 2e4],
        criteria,
    )
    @test valid_result.classification == :valid

    # Non-monotonic curve should be invalid
    invalid_result = classify_polarization_curve(
        [1.05, 0.92, 0.80, 0.85, 0.50],
        [0.0, 5e3, 1e4, 1.5e4, 2e4],
        criteria,
    )
    @test invalid_result.classification == :invalid

    # Full pipeline with IRD skipped
    result = run_validity_analysis(cfg)
    @test result isa ValidityAnalysisResult
    @test result.validation_summary.total_simulations == cfg.n_samples
    @test haskey(result.output_files, :parameter_classification_csv)
    @test haskey(result.output_files, :bounds_initial_yaml)
end

@testset "Parametrisation.SobolSensitivityAnalysis" begin
    using AlphaPEM.Parametrisation.SobolSensitivityAnalysis

    # Region construction
    regions = build_regions((0.4, 1.6))
    @test length(regions) == 3
    @test regions[1].name == :activation
    @test regions[2].name == :ohmic
    @test regions[3].name == :mass_transport

    # Regional AUC
    ifc = [0.0, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5]
    Ucell = [1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.1]
    aucs = compute_regional_aucs(ifc, Ucell, regions)
    @test length(aucs) == 3
    @test all(aucs .>= 0.0)

    # Input parameter construction
    cfg = SobolAnalysisConfig(
        fuel_cell_type = :ZSW_nominal,
        voltage_zone = :before_voltage_drop,
        N = 4,
        include_operating_conditions = true,
    )
    params = build_input_parameters(cfg)
    @test length(params) > 0
    @test all(p -> p.source in (:physical, :operating), params)

    # Excluded operating conditions
    cfg_excl = SobolAnalysisConfig(
        fuel_cell_type = :ZSW_nominal,
        voltage_zone = :before_voltage_drop,
        N = 4,
        include_operating_conditions = true,
        excluded_operating_conditions = [:Sa, :Sc],
    )
    params_excl = build_input_parameters(cfg_excl)
    @test :Sa ∉ [p.name for p in params_excl]
    @test :Sc ∉ [p.name for p in params_excl]
    @test :T_des ∈ [p.name for p in params_excl]

    # Sobol design matrices
    A, B = generate_sobol_design_matrices(cfg, params)
    @test size(A, 1) == length(params)
    @test size(A, 2) == cfg.N
    @test size(B) == size(A)

    # Operating condition validation
    oc = AlphaPEM.Config.OperatingConditions()
    @test is_valid_operating_conditions(oc)

    # Tiny end-to-end Sobol analysis (very few samples, no second order)
    cfg_tiny = SobolAnalysisConfig(
        fuel_cell_type = :ZSW_nominal,
        year = 2024,
        voltage_zone = :before_voltage_drop,
        N = 2,
        second_order = false,
        include_operating_conditions = false,
        nb_gc = 1,
        parallel = false,
        max_run_time_s = 30.0,
        output_dir = mktempdir(),
        save_curves = false,
    )
    result = run_sobol_analysis(cfg_tiny)
    @test result isa SobolAnalysisResult
    @test length(result.regions) == 3
    @test haskey(result.sobol_indices, :activation)
    @test haskey(result.sobol_indices, :ohmic)
    @test haskey(result.sobol_indices, :mass_transport)
end

@testset "Parametrisation.Calibration" begin
    using AlphaPEM.Parametrisation.Calibration
    using AlphaPEM.Config: SimulationConfig, PolarizationCalibrationParams, GAConfig, CalibrationConfig

    # GA and calibration configs
    ga_cfg = GAConfig(
        num_generations = 2,
        pop_size = 4,
        target_error = 5/100,
    )
    @test ga_cfg.pop_size == 4

    sim_cfg = SimulationConfig(
        type_fuel_cell = :ZSW_nominal,
        year = 2024,
        type_current = PolarizationCalibrationParams(),
        numerical_parameters = quick_numerical_params(),
        type_display = :no_display,
    )

    calib_cfg = CalibrationConfig(
        simulation_configs = [sim_cfg],
        ga_config = ga_cfg,
        parallel = false,
        output_dir = mktempdir(),
    )
    @test length(calib_cfg.simulation_configs) == 1

    # Parameter bounds for the fuel cell type
    using AlphaPEM.Parametrisation: ParametrisationCommon
    pb = ParametrisationCommon.bounds_for_fuel_cell(:ZSW_nominal, :full)
    @test !isempty(pb.bounds)
    @test all(b -> b.min < b.max, pb.bounds)
end

@testset "Parametrisation.CVExtraction" begin
    using AlphaPEM.Parametrisation.CVExtraction

    cv_file = joinpath(dirname(@__DIR__), "data", "ZSW", "2026", "cv", "cell_10.txt")

    cfg = CVExtractionConfig(
        area_cm2 = 283.87,
        cycle_for_extraction = 0,
        ignore_cycles_for_mean = [1],
        double_layer_limit_min = 0.30,
        double_layer_limit_max = 0.50,
        ohmic_drop_limit_min = 0.35,
        ohmic_drop_limit_max = 0.45,
        ecsa_limit_min = 0.09,
        ecsa_limit_max = 0.50,
        conversion_factor_uc_cm2 = 210.0,
        covering_degree = 0.77,
    )

    result = extract_cv_parameters(cv_file, cfg)

    @test result.cv_file_name == "cell_10"
    @test result.ecsa_adsorption_cm2 > 0.0
    @test result.ecsa_desorption_cm2 > 0.0
    @test result.crossover_a_cm2 > 0.0
    @test result.dlc_f_cm2 > 0.0
    @test result.scan_rate_vs > 0.0
    @test length(result.mean_cycle.U) > 0
    @test length(result.raw_mean_cycle.U) > 0
    @test length(result.anodic.U) > 0
    @test length(result.cathodic.U) > 0

    # Reference case: the values below are those reported by the original ZSW
    # Matlab tool for this file with the settings above and the 280 cm² area from its header.
    ref_file = joinpath(dirname(@__DIR__), "data", "reference_cv_for_checking_results.txt")
    ref_cfg = CVExtractionConfig(
        area_cm2 = 280.0,
        cycle_for_extraction = 0,
        ignore_cycles_for_mean = [1],
        double_layer_limit_min = 0.30,
        double_layer_limit_max = 0.50,
        ohmic_drop_limit_min = 0.35,
        ohmic_drop_limit_max = 0.45,
        ecsa_limit_min = 0.09,
        ecsa_limit_max = 0.50,
        conversion_factor_uc_cm2 = 210.0,
        covering_degree = 0.77,
    )

    ref = extract_cv_parameters(ref_file, ref_cfg)

    @test isapprox(ref.scan_rate_vs, 0.019910; rtol = 1e-4)   # 20 mV/s
    @test isapprox(ref.ecsa_desorption_cm2, 75.0; rtol = 0.01)
    @test isapprox(ref.ecsa_adsorption_cm2, 117.0; rtol = 0.01)
    @test isapprox(ref.crossover_a_cm2, 3.407e-3; rtol = 0.01)
    @test isapprox(ref.dlc_f_cm2, 4.6e-3; rtol = 0.02)
end

@testset "Utils" begin
    using AlphaPEM.Utils

    # Harmonic mean
    @test hmean(2.0, 2.0) ≈ 2.0
    @test hmean([2.0, 4.0]) ≈ 8.0 / 3.0

    # Physical property helpers
    @test rho_H2O_l(300.0) > 900.0
    @test Psat(300.0) > 0.0
    @test C_v_sat(300.0) > 0.0
    @test nu_l(300.0) > 0.0
end
