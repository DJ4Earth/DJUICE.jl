using dJUICE
using Test


@testset "Testing dJUICE" begin
# include("cost.jl")
# include("test.jl")
# include("test101.jl")
# include("test201.jl")
# include("test208.jl")
# include("test301.jl")
# include("test501.jl")
# include("testad.jl")
end

@testset "Model sturct tests" begin
	md = model()
	@test (typeof(md.mesh)) <: dJUICE.AbstractMesh
	@test (typeof(md.mesh)) == dJUICE.Mesh2dTriangle
	@test (typeof(md.friction)) <: dJUICE.AbstractFriction
	@test (typeof(md.friction)) <: dJUICE.BuddFriction

	md = model(md; friction=DNNFriction())
	@test (typeof(md.friction)) == dJUICE.DNNFriction
end
