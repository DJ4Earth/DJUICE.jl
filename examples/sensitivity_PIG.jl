using DJUICE
using MAT
using GLMakie

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, ".","PIG_Control_drag_dJUICE.mat"))

mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster 
md.stressbalance.maxiter = 100

#Now call AD!
md.inversion.iscontrol = 1
md.friction.coefficient = 30*ones(size(md.friction.coefficient))
md.inversion.independent = md.friction.coefficient
md.inversion.independent_string = "FrictionCoefficient"
md.inversion.dependent_string = ["SurfaceAbsVelMisfit"]

md = solve(md, :sb)

# the gradient
g = md.results["StressbalanceSolution"]["Gradient"]

plotmodel(md, abs.(g), caxis=(0.0,5e-3))

# save results to MAT for postproc
filename = "sensitivity_PIG.mat"
matwrite(filename, Dict("Gradient1" => g))
