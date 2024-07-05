using DJUICE
include("utils.jl")

md = model()
md = triangle(md,issmdir()*"/test/Exp/Pig.exp",20000.)
md = setmask( md,issmdir()*"/test/Exp/PigShelves.exp",issmdir()*"/test/Exp/PigIslands.exp")

#Initial velocity and geometry
x         = archread(issmdir()*"/test/Data/Pig.arch","x")
y         = archread(issmdir()*"/test/Data/Pig.arch","y")
vx_obs    = archread(issmdir()*"/test/Data/Pig.arch","vx_obs")
vy_obs    = archread(issmdir()*"/test/Data/Pig.arch","vy_obs")
index     = Int.(archread(issmdir()*"/test/Data/Pig.arch","index"))
surface   = archread(issmdir()*"/test/Data/Pig.arch","surface")
thickness = archread(issmdir()*"/test/Data/Pig.arch","thickness")
bed       = archread(issmdir()*"/test/Data/Pig.arch","bed")
md.initialization.vx=InterpFromMeshToMesh2d(index, x, y, vx_obs, md.mesh.x, md.mesh.y, 0.0)
md.initialization.vy=InterpFromMeshToMesh2d(index, x, y, vy_obs, md.mesh.x, md.mesh.y, 0.0)
md.geometry.surface = InterpFromMeshToMesh2d(index, x, y, surface, md.mesh.x, md.mesh.y, 0.0)
md.geometry.thickness = InterpFromMeshToMesh2d(index, x, y, thickness, md.mesh.x, md.mesh.y, 0.0)
md.geometry.base=md.geometry.surface .- md.geometry.thickness
md.geometry.bed =md.geometry.base
pos = findall(md.mask.ocean_levelset.<0)
md.geometry.bed[pos] = InterpFromMeshToMesh2d(index, x, y, bed, md.mesh.x[pos], md.mesh.y[pos])

md.materials.rheology_B=1.815730284801701e+08*ones(md.mesh.numberofvertices)
md.materials.rheology_n=3*ones(md.mesh.numberofelements)
md.friction.coefficient=50*ones(md.mesh.numberofvertices)
md.friction.p=ones(md.mesh.numberofvertices)
md.friction.q=ones(md.mesh.numberofvertices)

md.stressbalance.restol=0.05
md.stressbalance.reltol=1.0
md.stressbalance.abstol=NaN

#Boundary conditions
pos = findall(vec(sum(Int64.(md.mask.ocean_levelset[md.mesh.elements].<0), dims=2)) .> 0.0)
vertexonfloatingice=zeros(md.mesh.numberofvertices)
vertexonfloatingice[md.mesh.elements[pos,:]] .= 1
nodefront=(md.mesh.vertexonboundary .& (vertexonfloatingice.>0))
md.mask.ice_levelset[findall(nodefront)] .= 0

md.stressbalance.spcvx = NaN*ones(md.mesh.numberofvertices)
md.stressbalance.spcvy = NaN*ones(md.mesh.numberofvertices)
segmentsfront=md.mask.ice_levelset[md.mesh.segments[:,1:2]]==0
segments = findall(vec(sum(Int64.(md.mask.ice_levelset[md.mesh.segments[:,1:2]].==0), dims=2)) .!=2)
pos=md.mesh.segments[segments,1:2]
md.stressbalance.spcvx[pos] .= 0.0
md.stressbalance.spcvy[pos] .= 0.0

md=solve(md, :Stressbalance)

field_names =["Vx","Vy","Vel"]
field_tolerances=[NaN,NaN,NaN]
field_values= [(md.results["StressbalanceSolution"]["Vx"]),
               (md.results["StressbalanceSolution"]["Vy"]),
               (md.results["StressbalanceSolution"]["Vel"]) ]
compareArchive(@__FILE__, field_names, field_tolerances, field_values, :test)
