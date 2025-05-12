using DJUICE
include("utils.jl")

md = model()
md = triangle(md,issmdir()*"/test/Exp/Square.exp", 200000.)
md = setmask(md,"all","")

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
md.geometry.bed       = md.geometry.base .-500

#Initial velocity
x     = archread(issmdir()*"/test/Data/SquareShelf.arch","x")
y     = archread(issmdir()*"/test/Data/SquareShelf.arch","y")
vx    = archread(issmdir()*"/test/Data/SquareShelf.arch","vx")
vy    = archread(issmdir()*"/test/Data/SquareShelf.arch","vy")
index = Int.(archread(issmdir()*"/test/Data/SquareShelf.arch","index"))
md.initialization.vx=InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y,0.0)
md.initialization.vy=InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y,0.0)

md.materials.rheology_B=1.815730284801701e+08*ones(md.mesh.numberofvertices)
md.materials.rheology_n=3*ones(md.mesh.numberofelements)
md.friction.coefficient=20*ones(md.mesh.numberofvertices)
md.friction.p=ones(md.mesh.numberofvertices)
md.friction.q=ones(md.mesh.numberofvertices)

md.stressbalance.restol=0.10
md.stressbalance.reltol=0.02
md.stressbalance.abstol=NaN

#Boundary conditions
md.stressbalance.spcvx = NaN*ones(md.mesh.numberofvertices)
md.stressbalance.spcvy = NaN*ones(md.mesh.numberofvertices)
md.masstransport.spcthickness = NaN*ones(md.mesh.numberofvertices)
pos = findall(md.mesh.vertexonboundary)
md.stressbalance.spcvx[pos] .= 0.0
md.stressbalance.spcvy[pos] .= 0.0
md.masstransport.spcthickness[pos] .= md.geometry.thickness[pos]

md = SetIceShelfBC(md,issmdir()*"/test/Exp/SquareFront.exp")
# surface mass balance and basal melting
md.smb.mass_balance=10*ones(md.mesh.numberofvertices)
md.basalforcings.floatingice_melting_rate=5*ones(md.mesh.numberofvertices)
md.basalforcings.groundedice_melting_rate=5*ones(md.mesh.numberofvertices)

# mask
Lx = xmax - xmin
alpha = 2.0/3.0
md.mask.ice_levelset = ((md.mesh.x .- alpha*Lx).>0) .- ((md.mesh.x .- alpha*Lx).<0)

md.levelset.kill_icebergs = 0
# time stepping
md.timestepping.time_step = 10;
md.timestepping.final_time = 30;

# Transient
md.transient.isstressbalance=1;
md.transient.ismasstransport=1;
md.transient.issmb=1;
md.transient.ismovingfront=0;

md.calving.calvingrate=0*ones(md.mesh.numberofvertices)
md.frontalforcings.meltingrate=10000*ones(md.mesh.numberofvertices)
md.frontalforcings.ablationrate=10000*ones(md.mesh.numberofvertices)
md.levelset.spclevelset=NaN*ones(md.mesh.numberofvertices)
md.levelset.spclevelset[pos] = md.mask.ice_levelset[pos]

md=solve(md,:Transient)
field_names =["Vx1","Vy1","Vel1"]
field_tolerances=[4e-13,4e-13,4e-13]
field_values= [(md.results["TransientSolution"][1]["Vx"]),
					(md.results["TransientSolution"][1]["Vy"]),
					(md.results["TransientSolution"][1]["Vel"]) ]
#compareArchive(@__FILE__, field_names, field_tolerances, field_values, :test)

