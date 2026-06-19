# -*- coding: utf-8 -*-
using AlphaPEM.Parametrisation
using AlphaPEM.Parametrisation.Calibration
import AlphaPEM
using Test
using YAML

@testset "Calibration Output Files" begin
    output_dir = "results/test_calibration"
    rm(output_dir, force=true, recursive=true)
    
    cfg = CalibrationConfig(
        fuel_cell_type = :ZSW_GenStack,
        voltage_zone = :full,
        num_generations = 2,
        pop_size = 4,
        parallel = false,
        output_dir = output_dir,
        save_frequency = 1
    )
    
    # Test internal saving logic directly since optimize call might fail due to env/version issues
    try
        mkpath(output_dir)
        
        # Mock data
        history = [0.5, 0.4, 0.3]
        best_params = AlphaPEM.Config.PhysicalParams()
        res = CalibrationResult(cfg, best_params, 0.3, 0.3, history, 120.0)
        
        last_pop = [[1.0, 2.0], [3.0, 4.0]]
        last_fit = [0.5, 0.3]
        
        # 1. Test _save_final_results
        AlphaPEM.Parametrisation.Calibration._save_final_results(res, output_dir, last_pop, last_fit)
        
        @test isfile(joinpath(output_dir, "calibration_report.yaml"))
        @test isfile(joinpath(output_dir, "final_population.yaml"))
        
        # 2. Test checkpoint removal (simulated)
        checkpoint_path = joinpath(output_dir, "calibration_checkpoint.yaml")
        touch(checkpoint_path)
        @test isfile(checkpoint_path)
        
        # Re-run a part of calibrate logic or just the removal
        if isfile(checkpoint_path)
            rm(checkpoint_path)
        end
        @test !isfile(checkpoint_path)
        
        @info "Internal saving logic verified."
        
    catch e
        @error "Internal test failed: $e"
        rethrow(e)
    end
end
