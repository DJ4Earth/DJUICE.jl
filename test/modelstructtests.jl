module ModelStruct

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

@testset "Triangle" begin
	md = model()
	c = dJUICE.ExpStruct()
	c.name = "domainoutline"
	c.nods = 5
	c.density = 1.0
	c.x = [0.0, 1.0e6, 1.0e6, 0.0, 0.0]
	c.y = [0.0, 0.0, 1.0e6, 1.0e6, 0.0]
	c.closed = true
	contour = [c]
	md = triangle(md,contour,50000.) 
	@test md.mesh.numberofvertices == 340
	@test md.mesh.numberofelements == 614
end

end
