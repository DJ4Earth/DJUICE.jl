using Enzyme
using DJUICE
using MAT
using Test

using Optimization, OptimizationOptimJL

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMALL model (35 elements)
mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster 
md.stressbalance.maxiter = 20

#Now call AD!
md.inversion.iscontrol = 1
md.inversion.onlygrad = 1
md.inversion.independent = md.friction.coefficient
md.inversion.min_parameters = ones(md.mesh.numberofvertices)*(0.0)
md.inversion.max_parameters = ones(md.mesh.numberofvertices)*(1.0e3)
md.inversion.independent_string = "FrictionCoefficient"
md.inversion.dependent_string = ["SurfaceAbsVelMisfit"]

α = md.inversion.independent
∂J_∂α = zero(α)

femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
n = length(α)

DJUICE.costfunction(α, femmodel)
# test Enzyme autodiff only
dfemmodel = Enzyme.Compiler.make_zero(Base.Core.Typeof(femmodel), IdDict(), femmodel)
autodiff(set_runtime_activity(Enzyme.Reverse), DJUICE.costfunction, Active, Duplicated(α, ∂J_∂α), Duplicated(femmodel,dfemmodel))

# use user defined grad, errors!
#optprob = OptimizationFunction(DJUICE.costfunction, Optimization.AutoEnzyme())
#prob = Optimization.OptimizationProblem(optprob, α, femmodel, lb=md.inversion.min_parameters, ub=md.inversion.max_parameters)
#prob = Optimization.OptimizationProblem(optprob, α, femmodel)
#sol = Optimization.solve(prob,  Optimization.LBFGS())
#sol = Optimization.solve(prob, Optim.GradientDescent())
#sol = Optimization.solve(prob, Optim.NelderMead())
