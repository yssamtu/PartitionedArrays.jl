using Test

@testset "MatrixAssemblyBenchmarks" begin
    @testset "benchmark_tests" begin
        include("benchmark_tests.jl")
    end
end