using Enzyme
using DJUICE
using MAT
using Test

using ManualNLPModels
using MadNLP

onlygrad = 0

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "..", "data","PIG_testopt_small.mat")) #SMALL model (35 elements)
mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster
md.stressbalance.maxiter = 20


#Now call AD!
md.inversion.iscontrol = 1
md.inversion.onlygrad = onlygrad
md.inversion.independent = md.friction.coefficient
md.inversion.min_parameters = ones(md.mesh.numberofvertices)*(1.0)
md.inversion.max_parameters = ones(md.mesh.numberofvertices)*(200)
md.inversion.independent_string = "FrictionCoefficient"
md.inversion.dependent_string = ["SurfaceAbsVelMisfit"]

md.verbose.convergence = false

if onlygrad == 1
	md = solve(md, :sb)
end
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
	dx = zero(x)
	Enzyme.autodiff(set_runtime_activity(Enzyme.Reverse), DJUICE.costfunction, Active, Duplicated(x, dx), Duplicated(femmodel,dfemmodel))
	gx .= dx
end

nlp = NLPModel(
					α,
					f,
					grad = g!,
					lvar = md.inversion.min_parameters,
					uvar = md.inversion.max_parameters,
					)

#output = lbfgs(nlp)

results_qn = madnlp(
						  nlp;
						  linear_solver=LapackCPUSolver,
						  hessian_approximation=MadNLP.CompactLBFGS,
						  tol=1e-5,
						  max_iter=40,
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
