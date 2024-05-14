
using dJUICE
using MAT
using Test
using Enzyme
Enzyme.API.typeWarning!(false)
Enzyme.Compiler.RunAttributor[] = false

using Optimization

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

α = md.inversion.independent
α[2] = 10.0
∂J_∂α = zero(α)

femmodel=dJUICE.ModelProcessor(md, :StressbalanceSolution)
n = length(α)
# use user defined grad, errors!
optprob = OptimizationFunction(DJUICE.costfunction, Optimization.AutoEnzyme())
prob = Optimization.OptimizationProblem(optprob, α, femmodel, lb=md.inversion.min_parameters, ub=md.inversion.max_parameters, maxiters=1000)
#prob = Optimization.OptimizationProblem(optprob, α, femmodel)
#sol = Optimization.solve(prob, Optim.NelderMead())
sol = Optimization.solve(prob, Optimization.LBFGS())

# compute the gradient by enzyme
md = DJUICE.solve(md, :grad)
enz_grad = md.results["StressbalanceSolution"]["Gradient"]

# evaluating the gradient
sol.cache.f.grad(∂J_∂α, prob.u0, prob.p)

#independent_enum = DJUICE.StringToEnum(md.inversion.independent_string)
#DJUICE.InputUpdateFromVectorx(femmodel, sol.u, independent_enum, DJUICE.VertexSIdEnum)
#DJUICE.RequestedOutputsx(femmodel, [independent_enum])
#DJUICE.OutputResultsx(femmodel, md, :StressbalanceSolution)

#md = solve(md, :sb)

# compute gradient by finite differences at each node
#addJ = md.results["StressbalanceSolution"]["Gradient"]

