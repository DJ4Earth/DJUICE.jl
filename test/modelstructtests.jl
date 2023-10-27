using dJUICE
using Test
using MAT

@testset "Model basic sturct" begin
	md = model()
	@test (typeof(md.mesh)) <: dJUICE.AbstractMesh
	@test (typeof(md.mesh)) == dJUICE.Mesh2dTriangle
	@test (typeof(md.friction)) <: dJUICE.AbstractFriction
	@test (typeof(md.friction)) <: dJUICE.BuddFriction

	md = model(md; friction=DNNFriction())
	@test (typeof(md.friction)) == dJUICE.DNNFriction

	md = model(md; friction=SchoofFriction())
	@test (typeof(md.friction)) == dJUICE.SchoofFriction
end

@testset "Enums <-> String" begin
	@test dJUICE.StringToEnum("FrictionCoefficient") == dJUICE.FrictionCoefficientEnum
	@test dJUICE.StringToEnum("MaterialsRheologyB") == dJUICE.MaterialsRheologyBEnum
	@test dJUICE.EnumToString(dJUICE.FrictionCoefficientEnum) == "FrictionCoefficient"
	@test dJUICE.EnumToString(dJUICE.MaterialsRheologyBEnum) == "MaterialsRheologyB"
end

# quick test for AD
@testset "Cost function" begin
	file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMALL model (35 elements)
	mat  = read(file, "md")
	close(file)
	md = model(mat)

	#Now call AD!
	md.inversion.iscontrol = 1
	md.inversion.independent = md.friction.coefficient
	md.inversion.independent_string = "FrictionCoefficient"

	α = md.inversion.independent
	femmodel=dJUICE.ModelProcessor(md, :StressbalanceSolution)
	J1 = dJUICE.costfunction(femmodel, α)
	@test ~isnothing(J1)
end
