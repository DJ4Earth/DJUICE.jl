using DJUICE

md = model()
md = triangle(md,issmdir()*"/test/Exp/Square.exp",50000.)
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
md.geometry.bed       = md.geometry.base .-10

#Initial velocity
x     = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","x")
y     = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","y")
vx    = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","vx")
vy    = archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","vy")
index = Int.(archread(issmdir()*"/test/Data/SquareShelfConstrained.arch","index"))
md.initialization.vx=InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y,0.0)
md.initialization.vy=InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y,0.0)

md.materials.rheology_B=1.815730284801701e08*ones(md.mesh.numberofvertices)
md.materials.rheology_n=3*ones(md.mesh.numberofelements)
md.friction.coefficient=0.0*ones(md.mesh.numberofvertices)
md.friction.p=ones(md.mesh.numberofvertices)
md.friction.q=ones(md.mesh.numberofvertices)

md.stressbalance.restol=0.05
md.stressbalance.reltol=0.05
md.stressbalance.abstol=NaN

#Boundary conditions
md.stressbalance.spcvx = NaN*ones(md.mesh.numberofvertices)
md.stressbalance.spcvy = NaN*ones(md.mesh.numberofvertices)
pos = findall(md.mesh.vertexonboundary)
md.stressbalance.spcvx[pos] .= 0.0
md.stressbalance.spcvy[pos] .= 0.0

md=solve(md,:Stressbalance)

# Fields and tolerances to track changes
field_names =["Vx","Vy","Vel"]
field_tolerances=[4e-13,4e-13,4e-13]
field_values= [(md.results["StressbalanceSolution"]["Vx"]),
					(md.results["StressbalanceSolution"]["Vy"]),
					(md.results["StressbalanceSolution"]["Vel"]) ]
