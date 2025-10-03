using DJUICE
using MAT
using GLMakie
using JLD2

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "./Models/","PIG_Control_drag_dJUICE.mat"))

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

# save md
@save "./Models/PIG_Sensitivity.jld2" md

# the gradient
g = md.results["StressbalanceSolution"]["Gradient"]

# save results to MAT for postproc
filename = "./Models/PIG_djuice_gradient.mat"
matwrite(filename, Dict("Gradient1" => g))
