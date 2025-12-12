# DNN friction that learns y=x^1/3 and compare with Weertman friction
using DJUICE
using Lux
using JLD2

md = model()
md = model(md, friction=DNNFriction())
md = triangle(md,issmdir()*"/test/Exp/Square.exp",150000.)
md = setmask(md,"","")

#Geometry
hmin=300.
hmax=1000.
ymin=minimum(md.mesh.y)
ymax=maximum(md.mesh.y)
xmin=minimum(md.mesh.x)
xmax=maximum(md.mesh.x)
md.geometry.thickness = hmax .+ (hmin-hmax)*(md.mesh.y .- ymin)./(ymax-ymin) .+ 0.1*(hmin-hmax)*(md.mesh.x .- xmin)./(xmax-xmin)
md.geometry.base      = -md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness .+20
md.geometry.surface   = md.geometry.base+md.geometry.thickness
md.geometry.bed       = md.geometry.base

#Initial velocity
md.initialization.vx=zeros(md.mesh.numberofvertices)#InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y)
md.initialization.vy=zeros(md.mesh.numberofvertices)#InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y)

md.materials.rheology_B=1.815730284801701e+08*ones(md.mesh.numberofvertices)
md.materials.rheology_n=3*ones(md.mesh.numberofelements)

# load NN
@load issmdir()*"/test/Data/friction_DNN.jld2" trained_model
md.friction.model = trained_model.model.model
md.friction.ps = trained_model.ps
md.friction.st = trained_model.st
md.friction.C = 20*ones(md.mesh.numberofvertices)
md.friction.input_scale = 4e3/md.constants.yts
md.friction.output_scale = (md.friction.input_scale)^(1.0/3.0)

md.stressbalance.restol=0.05
md.stressbalance.reltol=0.05
md.stressbalance.abstol=NaN

#Boundary conditions
md.stressbalance.spcvx = NaN*ones(md.mesh.numberofvertices)
md.stressbalance.spcvy = NaN*ones(md.mesh.numberofvertices)
pos = findall(md.mesh.vertexonboundary)
md.stressbalance.spcvx[pos] .= 0.0
md.stressbalance.spcvy[pos] .= 0.0


# Weertman
md2 = model(deepcopy(md), friction=WeertmanFriction())
md2.friction.C = 20*ones(md.mesh.numberofvertices)
md2.friction.m = 3*ones(md.mesh.numberofvertices)

md = solve(md,:Stressbalance)
md2 = solve(md2,:Stressbalance)

@testset "DNNFriction for Weertman" begin
	err = abs.(md2.results["StressbalanceSolution"]["Vel"]-md.results["StressbalanceSolution"]["Vel"])./md.results["StressbalanceSolution"]["Vel"]
	@test sum(filter(!isnan, err)) < 1e-2
end
