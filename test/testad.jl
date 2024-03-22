module enzymeDiff_grad_frictionC

using dJUICE
using MAT
using Test
using Enzyme

Enzyme.Compiler.RunAttributor[] = false

#Load model from MATLAB file
#file = matopen(joinpath(@__DIR__, "..", "data","temp12k.mat")) #BIG model
file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMALL model (35 elements)
mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster 
md.stressbalance.maxiter = 20

#Now call AD!
md.inversion.iscontrol = 1
md.inversion.independent = md.friction.coefficient
md.inversion.independent_string = "FrictionCoefficient"

md = solve(md, :grad)

# compute gradient by finite differences at each node
addJ = md.results["StressbalanceSolution"]["Gradient"]


@testset "Quick AD test with Cost function" begin
   #Now call AD!
   md.inversion.iscontrol = 1
   md.inversion.independent = md.friction.coefficient
   md.inversion.independent_string = "FrictionCoefficient"

   α = md.inversion.independent
   femmodel=dJUICE.ModelProcessor(md, :StressbalanceSolution)
   J1 = dJUICE.costfunction(femmodel, α)
   @test ~isnothing(J1)
end

@testset "AD gradient calculation for FrictionC" begin
	α = md.inversion.independent
	delta = 1e-7
	femmodel=dJUICE.ModelProcessor(md, :StressbalanceSolution)
	J1 = dJUICE.costfunction(femmodel, α)
	for i in 1:md.mesh.numberofvertices
		dα = zero(md.friction.coefficient)
		dα[i] = delta
		femmodel=dJUICE.ModelProcessor(md, :StressbalanceSolution)
		J2 = dJUICE.costfunction(femmodel, α+dα)
		dJ = (J2-J1)/delta

		@test abs(dJ - addJ[i])< 1e-5
	end
end

end
