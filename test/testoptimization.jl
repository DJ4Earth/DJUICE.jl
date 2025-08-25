using Enzyme
using DJUICE
using MAT
using Test

#using Optimization, OptimizationNLopt
using ManualNLPModels, JSOSolvers
using MadNLP


#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMALL model (35 elements)
mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster
md.stressbalance.maxiter = 20

#Now call AD!
md.inversion.iscontrol = 1
md.inversion.onlygrad = 0
md.inversion.independent = md.friction.coefficient
md.inversion.min_parameters = ones(md.mesh.numberofvertices)*(0.0)
md.inversion.max_parameters = ones(md.mesh.numberofvertices)*(1.0e3)
md.inversion.independent_string = "FrictionCoefficient"
md.inversion.dependent_string = ["SurfaceAbsVelMisfit"]

α = md.inversion.independent
∂J_∂α = zero(α)


# cost function f
f(x) = begin
	femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
	DJUICE.costfunction(x, femmodel)
end

# Enzyme gradient g
g!(gx, x) = begin
	femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
	dfemmodel = Enzyme.Compiler.make_zero(Base.Core.Typeof(femmodel), IdDict(), femmodel)
	Enzyme.autodiff(set_runtime_activity(Enzyme.Reverse), DJUICE.costfunction, Active, Duplicated(x, gx), Duplicated(femmodel,dfemmodel))
	gx
end

nlp = NLPModel(
  α,
  f,
  grad = g!,
)

#output = lbfgs(nlp)

results_qn = madnlp(
           nlp;
           linear_solver=LapackCPUSolver,
           hessian_approximation=MadNLP.CompactLBFGS,
       )



# use user defined grad, errors!
#optprob = OptimizationFunction(DJUICE.costfunction, Optimization.AutoEnzyme(; mode=set_runtime_activity(Enzyme.Reverse)))
#prob = Optimization.OptimizationProblem(optprob, α, femmodel, lb=md.inversion.min_parameters, ub=md.inversion.max_parameters)
#prob = Optimization.OptimizationProblem(optprob, α, femmodel)
#sol = Optimization.solve(prob,  Optimization.LBFGS(), maxiters=10)
#sol = Optimization.solve(prob, Optim.GradientDescent(), maxiters=10)
#sol = Optimization.solve(prob, Optim.LBFGS(), maxiters=10)
#sol = Optimization.solve(prob, Optim.NelderMead())
#sol = Optimization.solve(prob, Optimization.LBFGS(), maxiters = 100)

#using JuMP, Optim
#m = Model(Optim.Optimizer);
#set_optimizer_attribute(m, "method", LBFGS())
#@variable(m, α)
#@variable(m, femmodel)
#@objective(m, Min, DJUICE.costfunction(α, femmodel))
