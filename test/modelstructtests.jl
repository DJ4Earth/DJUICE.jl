module ModelStruct

using DJUICE
using Test
using MAT

@testset "Model basic sturct" begin
	md = model()
	@test (typeof(md.mesh)) <: DJUICE.AbstractMesh
	@test (typeof(md.mesh)) == DJUICE.Mesh2dTriangle
	@test (typeof(md.friction)) <: DJUICE.AbstractFriction
	@test (typeof(md.friction)) <: DJUICE.BuddFriction

	md = model(md; friction=DNNFriction())
	@test (typeof(md.friction)) == DJUICE.DNNFriction

	md = model(md; friction=SchoofFriction())
	@test (typeof(md.friction)) == DJUICE.SchoofFriction
end

@testset "Enums <-> String" begin
	@test DJUICE.StringToEnum("FrictionCoefficient") == DJUICE.FrictionCoefficientEnum
	@test DJUICE.StringToEnum("MaterialsRheologyB") == DJUICE.MaterialsRheologyBEnum
	@test DJUICE.EnumToString(DJUICE.FrictionCoefficientEnum) == "FrictionCoefficient"
	@test DJUICE.EnumToString(DJUICE.MaterialsRheologyBEnum) == "MaterialsRheologyB"
end

@testset "Triangle" begin
	md = model()
	c = DJUICE.ExpStruct()
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

@testset "IceVolume" begin
	file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMALL model (35 elements)
	mat  = read(file, "md")
	close(file)
	md = model(mat)

	#make model run faster
	md.stressbalance.maxiter = 20
	md.verbose.convergence = false

	#Now call AD!
	md.inversion.iscontrol = 1
	md.inversion.independent = md.friction.coefficient
	md.inversion.independent_string = "FrictionCoefficient"
	md.inversion.dependent_string = ["IceVolume"]
	α = md.inversion.independent
	femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
	J1 = DJUICE.CostFunction(α, femmodel)
	intH = DJUICE.IntegrateOverDomain(md, md.geometry.thickness)
	@test isapprox(J1, intH)
end

@testset "Verbose" begin
	md = model()
	DJUICE.SetVerbosityLevel(md)
	@test DJUICE.verbositylevel == 0
	for f in fieldnames(typeof(md.verbose))
		setfield!(md.verbose,f,true)
	end
	DJUICE.SetVerbosityLevel(md)
	@test DJUICE.VerboseMProcessor() == true
	@test DJUICE.VerboseModule() == true
	@test DJUICE.VerboseSolution() == true
	@test DJUICE.VerboseSolver() == true
	@test DJUICE.VerboseConvergence() == true
	@test DJUICE.VerboseControl() == true
	@test DJUICE.VerboseAutodiff() == true
	md.verbose.solver = false
	DJUICE.SetVerbosityLevel(md)
	@test DJUICE.VerboseSolver() == false
end

end
