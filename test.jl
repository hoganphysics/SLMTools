using SLMTools


L1 = (1:10, 1:20)
L2 = dualShiftLattice(L1)  # Creating a dual lattice for comparison
flambda = 1.0

dpField = dualPhase(L1)
isa(dpField, LF{ComplexPhase})
# @test dpField.L == dualShiftLattice(L1)