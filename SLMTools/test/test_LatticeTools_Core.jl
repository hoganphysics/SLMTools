using Test
using SLMTools.LatticeTools.LatticeCore

@testset "Lattice tests" begin
    @testset "Function accepting Lattice" begin
        function test_function(l::Lattice{2})
            return sum(length(r) for r in l)
        end
        l = (1:10, 1:20)
        @test test_function(l) == 30
    end
end
@testset "natlat Tests" begin
    @testset "1D lattice tests" begin
        # Test for n = 5
        test_lattice = natlat(5)
        expected_lattice = [-2, -1, 0, 1, 2] ./ sqrt(5)
        @test test_lattice == expected_lattice
    end

    @testset "ND lattice tests" begin
        # Test for a 2D lattice
        test_lattice = natlat((4, 4))
        expected_lattice = (-1.0:0.5:0.5, -1.0:0.5:0.5)
        @test test_lattice == expected_lattice
    end
end