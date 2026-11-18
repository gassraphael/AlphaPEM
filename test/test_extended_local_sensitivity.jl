using AlphaPEM
using Test

@testset "ExtendedLocalSensitivityAnalysis sensitivity (eq. 7)" begin
    ELSA = AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis
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
