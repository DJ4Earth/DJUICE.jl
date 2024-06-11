using DJUICE
using MAT

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "..", "Data","PIG_Control_drag_dJUICE.mat"))

mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster 
md.stressbalance.maxiter = 20

#Now call AD!
md.inversion.iscontrol = 1
md.inversion.independent = md.friction.coefficient
md.inversion.independent_string = "FrictionCoefficient"

md = solve(md, :grad)

# the gradient
g = md.results["StressbalanceSolution"]["Gradient"]
