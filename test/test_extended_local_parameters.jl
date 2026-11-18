using AlphaPEM
using Test

@testset "ExtendedLocalSensitivityAnalysis parameters" begin
    ELSA = AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis
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

    # Subset and errors
    sub, _ = ELSA.select_parameters(ELSA.ExtendedLocalSensitivityConfig(parameter_subset = [:Hmem, :Re]))
    @test Set(p.name for p in sub) == Set([:Hmem, :Re])
    @test_throws ArgumentError ELSA.select_parameters(ELSA.ExtendedLocalSensitivityConfig(parameter_subset = [:nope]))
    bad = AlphaPEM.Parametrisation.Calibration.ParameterBound(:x, 0.0, 10.0, :real, :log, "", "")
    @test_throws ArgumentError ELSA.to_unit(bad, 1.0)
end
