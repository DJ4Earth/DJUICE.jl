# AD stressbalance rheology
using Enzyme

using DJUICE
using MAT
using Test

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
md.inversion.independent = md.materials.rheology_B
md.inversion.independent_string = "MaterialsRheologyB"
md.inversion.dependent_string = ["DragCoefficientAbsGradient"]

md = solve(md, :sb)

# compute gradient by finite differences at each node
addJ = md.results["StressbalanceSolution"]["Gradient"]

@testset "Quick AD test with Cost function" begin
   #Now call AD!
   md.inversion.iscontrol = 1
	md.inversion.independent = md.materials.rheology_B
	md.inversion.independent_string = "MaterialsRheologyB"

   α = md.inversion.independent
   femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
   J1 = DJUICE.costfunction(α, femmodel)
   @test ~isnothing(J1)
end

@testset "AD results RheologyB" begin
	α = md.inversion.independent
	delta = 1e-8
	femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
	J1 = DJUICE.costfunction(α, femmodel)
	for i in 1:md.mesh.numberofvertices
		dα = zero(md.friction.coefficient)
		dα[i] = delta
		femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
		J2 = DJUICE.costfunction(α+dα, femmodel)
		dJ = (J2-J1)/delta

		@test abs(dJ - addJ[i])< 1e-6
	end
end
