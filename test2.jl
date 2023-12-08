using SLMTools

data = rand(ComplexF64, 10, 10)
lattice = (1:10, 1:10)
flambda = 1.5
lf = LatticeField{ComplexPhase}(data, lattice)
sub_lf = lf[1:5, 1:5]

@test typeof(sub_lf) == LatticeField{ComplexPhase,ComplexF64,2}
@test size(sub_lf.data) == (5, 5)
@test sub_lf.L == (1:5, 1:5)









