using dJUICE
using Test

@testset "Model sturct tests" begin
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
