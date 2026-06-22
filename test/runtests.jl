using AlphaPEM
using Test

@testset "AlphaPEM.jl" begin
    # Basic loading test
    @test true
    
    # Verify access to core modules and types
    @test isdefined(AlphaPEM, :Config)
    @test isdefined(AlphaPEM.Config, :SimulationConfig)
    @test isdefined(AlphaPEM, :Core)
    @test isdefined(AlphaPEM, :Application)
end
