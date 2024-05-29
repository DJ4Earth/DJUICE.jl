using Enzyme
Enzyme.API.looseTypeAnalysis!(false)
Enzyme.API.strictAliasing!(false)
Enzyme.API.typeWarning!(false)
Enzyme.Compiler.RunAttributor[] = false

using Optimization, OptimizationOptimJL


function Control_Core(md::model, femmodel::FemModel) #{{{
	# solve for optimization
	# TODO: just a first try, need to add all the features
	α = md.inversion.independent
	∂J_∂α = zero(α)
	n = length(α)
	# use user defined grad, errors!
	#optprob = OptimizationFunction(costfunction, Optimization.AutoEnzyme(), grad=computeGradient(∂J_∂α, α, femmodel))
	optprob = OptimizationFunction(costfunction, Optimization.AutoEnzyme())
	prob = Optimization.OptimizationProblem(optprob, α, femmodel, lb=md.inversion.min_parameters, ub=md.inversion.max_parameters)
	sol = Optimization.solve(prob, Optim.LBFGS())

	independent_enum = StringToEnum(md.inversion.independent_string)
	InputUpdateFromVectorx(femmodel, sol.u, independent_enum, VertexSIdEnum)
	RequestedOutputsx(femmodel, [independent_enum])
end#}}}
function computeGradient(md::model, femmodel::FemModel) #{{{
	#independent variable
	α = md.inversion.independent
	#initialize derivative as 0
	∂J_∂α = zero(α)
	# Compute Gradient
	computeGradient(∂J_∂α, α, femmodel)

	#Put gradient in results
	InputUpdateFromVectorx(femmodel, ∂J_∂α, GradientEnum, VertexSIdEnum)
	RequestedOutputsx(femmodel, [GradientEnum])
end#}}}
function computeGradient(∂J_∂α::Vector{Float64}, α::Vector{Float64}, femmodel::FemModel) #{{{
	# zero ALL depth of the model, make sure we get correct gradient
	dfemmodel = Enzyme.Compiler.make_zero(Base.Core.Typeof(femmodel), IdDict(), femmodel)
	# compute the gradient
	autodiff(Enzyme.Reverse, costfunction, Duplicated(α, ∂J_∂α), Duplicated(femmodel,dfemmodel))
end#}}}
