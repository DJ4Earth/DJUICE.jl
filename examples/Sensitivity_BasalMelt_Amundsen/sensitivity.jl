using DJUICE
using MAT
using JLD2

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "./Models/", "Amundsen_Controls_dJUICE.mat"))

mat  = read(file, "md")
close(file)
md = model(mat, friction=DJUICE.SchoofFriction(), basalforcings=DJUICE.LinearBasalforcings())

md.timestepping.time_step = 0.1
md.timestepping.start_time = 1995
md.timestepping.final_time = 1995.5
#make model run faster 
#md.stressbalance.maxiter = 100

# manually set basal forcings
md.basalforcings.deepwater_melting_rate = 50.0
md.basalforcings.upperwater_melting_rate = 0.0
md.basalforcings.deepwater_elevation = -500.0
md.basalforcings.upperwater_elevation = 0.0
md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices)
md.basalforcings.perturbation_melting_rate = zeros(md.mesh.numberofvertices)

##Now call AD!
md.inversion.iscontrol = 1
md.inversion.independent = md.basalforcings.perturbation_melting_rate
md.inversion.independent_string = "BasalforcingsPerturbationMeltingRate"
md.inversion.dependent_string = ["IceVolumeAboveFloatation"]
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy

# solve
md = solve(md, :Transient)

# save
@save "./Models/Amundsen_Sensitivity.jld2" md

# the gradient
g = md.results["TransientSolution"][5]["Gradient"]

# save results to MAT for postproc
filename = "./Models/Amudsen_djuice_gradient.mat"
matwrite(filename, Dict("Gradient1" => g))
