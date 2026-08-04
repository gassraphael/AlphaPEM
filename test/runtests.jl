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
    @test cfg.type_fuel_cell == :ZSW_GenStack
    @test validate_config(cfg) === cfg

    # Custom configuration
    cfg_custom = SimulationConfig(
        type_fuel_cell = :ZSW_GenStack,
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
    fc = create_fuelcell(:ZSW_GenStack, :full)
    @test fc isa AbstractFuelCell
    @test fc.physical_parameters isa AlphaPEM.Config.PhysicalParams
    @test fc.operating_conditions isa AlphaPEM.Config.OperatingConditions

    # Other supported fuel cell types
    @test create_fuelcell(:EH31_2022, :full) isa AbstractFuelCell
    @test create_fuelcell(:ZSW_GenStack_T_84, :before_voltage_drop) isa AbstractFuelCell

    # Accessor functions
    @test physical_parameters(fc) isa AlphaPEM.Config.PhysicalParams
    @test operating_conditions(fc, :ZSW_GenStack) isa AlphaPEM.Config.OperatingConditions
    @test undetermined_parameters(fc, :full) isa Vector{Tuple{Symbol, Float64, Float64}}
    @test !isempty(undetermined_parameters(fc, :full))
end

@testset "Currents" begin
    using AlphaPEM.Currents
    using AlphaPEM.Fuelcell: create_fuelcell

    fc = create_fuelcell(:ZSW_GenStack, :full)

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
        type_fuel_cell = :ZSW_GenStack,
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
        type_fuel_cell = :ZSW_GenStack,
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
        type_fuel_cell = :ZSW_GenStack,
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
        type_fuel_cell = :ZSW_GenStack,
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
            type_fuel_cell = :ZSW_GenStack,
            type_current = pola_params,
            numerical_parameters = num_params,
            voltage_zone = :full,
            type_display = :no_display,
            display_timing = :postrun,
        ),
        SimulationConfig(
            type_fuel_cell = :ZSW_GenStack_T_84,
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
        type_fuel_cell = :ZSW_GenStack,
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
    pb = ParametrisationCommon.bounds_for_fuel_cell(:ZSW_GenStack, :full)
    @test !isempty(pb.bounds)
    @test all(b -> b.min < b.max, pb.bounds)
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
