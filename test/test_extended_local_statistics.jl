using AlphaPEM
using Test

@testset "ExtendedLocalSensitivityAnalysis statistics (eq. 8-9)" begin
    ELSA = AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

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
    ELSA = AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

    # Rows: 1 and 2 are collinear (same direction, different norm), 3 orthogonal to both, 4 is null
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
    @test ELSA.collinearity_index([1.0 -2.0; -2.0 4.0])[1, 2] ≈ 1.0   # opposite signs are still collinear
    @test 0.0 <= ELSA.collinearity_index([1.0 2.0 3.0; 2.0 -1.0 0.5])[1, 2] <= 1.0

    # Per point and median over points
    S_n = zeros(3, 2, 2)
    S_n[1, :, :] = [1.0 0.0; 1.0 0.0]               # ψ = 1
    S_n[2, :, :] = [1.0 0.0; 0.0 1.0]               # ψ = 0
    S_n[3, :, :] = [1.0 0.0; 1.0 0.0]               # ψ = 1
    Ψs = ELSA.collinearity_indices(S_n)
    @test size(Ψs) == (3, 2, 2)
    @test Ψs[:, 1, 2] ≈ [1.0, 0.0, 1.0]
    @test ELSA.median_collinearity(S_n)[1, 2] ≈ 1.0

    # Failed point (NaN) ignored
    S_nan = copy(S_n)
    S_nan[1, 1, :] .= NaN
    @test isnan(ELSA.collinearity_indices(S_nan)[1, 1, 2])
    @test ELSA.median_collinearity(S_nan)[1, 2] ≈ 0.5

    @test ELSA.influential_mask([0.5, 0.05, NaN, 0.1], 0.1) == BitVector([true, false, false, true])
    @test ELSA.influential_mask([Inf, -Inf, 0.2], 0.1) == BitVector([false, false, true])
    @test_throws ArgumentError ELSA.influential_mask([0.2], -0.1)
    @test_throws ArgumentError ELSA.influential_mask([0.2], Inf)
end
