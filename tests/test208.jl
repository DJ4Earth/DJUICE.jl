#!/Applications/Julia-1.6.app/Contents/Resources/julia/bin/julia
include("../issm.jl")
using .ISSM

md = model()
md = triangle(md,issmdir()*"/test/Exp/Square.exp", 150000.)
md = setmask(md, "all", "")

#Geometry
hmin=300.
hmax=1000.
ymin=minimum(md.mesh.y)
ymax=maximum(md.mesh.y)
xmin=minimum(md.mesh.x)
xmax=maximum(md.mesh.x)
md.geometry.thickness = hmax .+ (hmin-hmax)*(md.mesh.y .- ymin)./(ymax-ymin) .+ 0.1*(hmin-hmax)*(md.mesh.x .- xmin)./(xmax-xmin)
md.geometry.base      = -md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness
md.geometry.surface   = md.geometry.base+md.geometry.thickness
md.geometry.bed       = md.geometry.base .-10

#Initial velocity
x     = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","x")
y     = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","y")
vx    = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","vx")
vy    = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","vy")
index = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","index")
md.initialization.vx=zeros(md.mesh.numberofvertices)#InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y)
md.initialization.vy=zeros(md.mesh.numberofvertices)#InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y)

md.materials.rheology_B=1.815730284801701e+08*ones(md.mesh.numberofvertices)
md.materials.rheology_n=3*ones(md.mesh.numberofelements)
md.friction.coefficient=20*ones(md.mesh.numberofvertices)

md.stressbalance.restol=0.1
md.stressbalance.reltol=0.02
md.stressbalance.abstol=NaN
md.timestepping.start_time = 0.0
md.timestepping.final_time = 3.0
md.timestepping.time_step  = 1.0

#Boundary conditions
nodefront=ContourToNodes(md.mesh.x,md.mesh.y,issmdir()*"/test/Exp/SquareFront.exp",2.0) .& md.mesh.vertexonboundary
md.stressbalance.spcvx = NaN*ones(md.mesh.numberofvertices)
md.stressbalance.spcvy = NaN*ones(md.mesh.numberofvertices)
pos = findall(md.mesh.vertexonboundary .& .~nodefront)
md.mask.ice_levelset[findall(nodefront)] .= 0

segmentsfront=md.mask.ice_levelset[md.mesh.segments[:,1:2]]==0
segments = findall(vec(sum(Int64.(md.mask.ice_levelset[md.mesh.segments[:,1:2]].==0), dims=2)) .!=2)
pos=md.mesh.segments[segments,1:2]
md.stressbalance.spcvx[pos] .= 0.0
md.stressbalance.spcvy[pos] .= 0.0

md.smb.mass_balance = zeros(md.mesh.numberofvertices)
md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices)
md.basalforcings.floatingice_melting_rate = ones(md.mesh.numberofvertices)
md.masstransport.spcthickness = NaN*ones(md.mesh.numberofvertices)

md=solve(md,"Transient")
