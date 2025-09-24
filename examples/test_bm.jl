using DJUICE
using MAT
using JLD2
#using GLMakie

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, ".", "","ModelBig.mat"))

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

#
##Now call AD!
md.inversion.iscontrol = 1
#md.friction.coefficient = 30*ones(size(md.friction.coefficient))
#md.inversion.independent = md.friction.C
#md.inversion.independent_string = "FrictionC"

md.inversion.independent = md.basalforcings.perturbation_melting_rate
md.inversion.independent_string = "BasalforcingsPerturbationMeltingRate"
md.inversion.dependent_string = ["IceVolumeAboveFloatation"]
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy


#md = solve(md, :sb)
md = solve(md, :Transient)

# save
@save "test_bm_pm.jld2" md

# the gradient
#g = md.results["StressbalanceSolution"]["Gradient"]

#plotmodel(md, abs.(g), caxis=(0.0,5e-3))
