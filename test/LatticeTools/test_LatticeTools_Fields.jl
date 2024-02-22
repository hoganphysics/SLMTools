using Test
using SLMTools.LatticeTools

# Sample data and lattices for testing
data = rand(ComplexF64, 10, 10)
lattice = (1:10, 1:10)
flambda = 1.5

# Test LatticeField constructor
@testset "LatticeField Constructor Tests" begin
    lf = LatticeField{ComplexPhase}(data, lattice, flambda)
    @test size(lf.data) == (10, 10)
    @test lf.L == lattice
    @test lf.flambda ≈ 1.5

    # Test for DimensionMismatch
    data_mismatch = rand(ComplexF64, 8, 8)
    @test_throws DimensionMismatch LatticeField{ComplexPhase}(data_mismatch, lattice)
end

# Test convenience constructor
@testset "Convenience Constructor Tests" begin
    lf = LatticeField{ComplexPhase}(data, lattice)
    @test size(lf.data) == (10, 10)
    @test lf.L == lattice
    @test lf.flambda ≈ 1.0
end

# Test equality checking for lattices
@testset "Lattice Equality Tests" begin
    lattice2 = (1:10, 1:10)
    @test elq(lattice, lattice2) === nothing

    lattice3 = (1:9, 1:10)
    @test_throws DimensionMismatch elq(lattice, lattice3)

    lattice4 = (10, 10)|>natlat
    @test_throws DomainError elq(lattice, lattice4)
end

# Test minimal methods for LatticeField as array
@testset "LatticeField Array Behavior Tests" begin
    lf = LatticeField{ComplexPhase}(data, lattice)

    @test size(lf) == (10, 10)
    @test lf[1, 1] == data[1, 1]
    @test lf[5, 5] != 0

    # Modify an element
    lf[2, 2] = 0
    @test lf[2, 2] == 0
end

# Test subfield creation
@testset "Subfield Creation Tests" begin
    lf = LatticeField{ComplexPhase}(data, lattice)
    sub_lf = lf[1:5, 1:5]

    @test typeof(sub_lf) == LatticeField{ComplexPhase,ComplexF64,2}
    @test size(sub_lf.data) == (5, 5)
    @test sub_lf.L == (1:5, 1:5)
end

# Add more tests as needed
