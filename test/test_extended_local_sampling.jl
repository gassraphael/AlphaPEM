using AlphaPEM
using Test

@testset "ExtendedLocalSensitivityAnalysis sampling and screening (step 4)" begin
    ELSA = AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

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
        # candidate k has a constant bias θ[1] on every dataset; negative → failure
        phys = [0.1 0.5 -1.0 0.2; 0.0 0.0 0.0 0.0]
        sim(θ, d) = θ[1] < 0 ? nothing : (fill(0.7 + θ[1], 3), fill(0.7, 3))
        for par in (false, true)
            cfg = ELSA.ExtendedLocalSensitivityConfig(parallel = par, screening_datasets = ds,
                                                      screening_threshold_V = 0.3)
            rec = ELSA.screen_candidates(phys, ds, cfg; simulate_fn = sim)
            @test [r.status for r in rec] == [:accepted, :rejected, :failed, :accepted]
            @test rec[1].error_V ≈ 0.1
            @test length(rec[1].dataset_errors_V) == 2
            @test [r.id for r in rec] == 1:4
        end
        cfg = ELSA.ExtendedLocalSensitivityConfig(screening_datasets = ds)
        @test_throws ErrorException ELSA.screen_candidates(phys[:, 2:3], ds, cfg; simulate_fn = sim)

        # 10% exactly is allowed; anything below it aborts the screening.
        threshold_design = zeros(1, 10)
        pass_first(θ, d) = θ[1] == 1.0 ? (fill(0.7, 3), fill(0.7, 3)) :
                                         (fill(1.1, 3), fill(0.7, 3))
        threshold_design[1, 1] = 1.0
        @test length(ELSA.screen_candidates(threshold_design, ds, cfg; simulate_fn = pass_first)) == 10
        below_threshold_design = hcat(threshold_design, zeros(1, 1))
        @test_throws ErrorException ELSA.screen_candidates(
            below_threshold_design, ds, cfg; simulate_fn = pass_first)
    end

    @testset "config" begin
        cfg = ELSA.ExtendedLocalSensitivityConfig()
        @test cfg.n_candidates == 2^9
        @test cfg.include_nominal == false
        @test !hasproperty(cfg, :scramble) && !hasproperty(cfg, :n_target)
        @test ELSA.validate_config(ELSA.ExtendedLocalSensitivityConfig(n_candidates = 350)).n_candidates == 350
        @test_throws Exception ELSA.validate_config(ELSA.ExtendedLocalSensitivityConfig(screening_datasets = ELSA.ScreeningDataset[]))
    end
end
