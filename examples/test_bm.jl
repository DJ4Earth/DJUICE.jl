using DJUICE
using MAT
using GLMakie

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, ".", "","ModelBig.mat"))

mat  = read(file, "md")
close(file)
md = model(mat, friction=DJUICE.SchoofFriction())

md.timestepping.time_step = 0.1
#make model run faster 
#md.stressbalance.maxiter = 100
#
##Now call AD!
md.inversion.iscontrol = 1
#md.friction.coefficient = 30*ones(size(md.friction.coefficient))
md.inversion.independent = md.friction.C
md.inversion.independent_string = "FrictionC"
md.inversion.dependent_string = ["IceVolumeAboveFloatation"]
md.inversion.vx_obs = md.initialization.vx
md.inversion.vy_obs = md.initialization.vy


#md = solve(md, :sb)
md = solve(md, :Transient)

# the gradient
#g = md.results["StressbalanceSolution"]["Gradient"]

#plotmodel(md, abs.(g), caxis=(0.0,5e-3))
