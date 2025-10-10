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
md.inversion.maxiter = 40
md.inversion.tol = 1e-5

md.verbose.convergence = false

#md = solve(md, :sb);

α = md.inversion.independent
∂J_∂α = zero(α)

# cost function f
f(x) = begin
	femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
	DJUICE.CostFunction(x, femmodel)
end

# Enzyme gradient g
g!(gx, x) = begin
	femmodel=DJUICE.ModelProcessor(md, :StressbalanceSolution)
	DJUICE.ComputeGradient!(gx, x, femmodel)
end

nlp = NLPModel(
					md.inversion.independent,
					f,
					grad = g!,
					lvar = md.inversion.min_parameters,
					uvar = md.inversion.max_parameters,
					)

results_qn = madnlp(
						  nlp;
						  linear_solver=LapackCPUSolver,
						  hessian_approximation=MadNLP.CompactLBFGS,
						  tol=md.inversion.tol,
						  max_iter=md.inversion.maxiter,
						  )


