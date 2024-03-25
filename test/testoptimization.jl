module enzymeDiff

using dJUICE
using MAT
using Test
using Enzyme

Enzyme.Compiler.RunAttributor[] = false

#Load model from MATLAB file
#file = matopen(joinpath(@__DIR__, "..", "data","temp12k.mat")) #BIG model
file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMALL model (35 elements)
mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster 
md.stressbalance.maxiter = 20

#Now call AD!
md.inversion.iscontrol = 1
md.inversion.independent = md.friction.coefficient
md.inversion.min_parameters = ones(md.mesh.numberofvertices)*(0.0)
md.inversion.max_parameters = ones(md.mesh.numberofvertices)*(1.0e3)
md.inversion.independent_string = "FrictionCoefficient"

md = solve(md, :sb)

# compute gradient by finite differences at each node
#addJ = md.results["StressbalanceSolution"]["Gradient"]

end
