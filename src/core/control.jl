using Enzyme
using ManualNLPModels
using MadNLP

function Control_Core(md::model, femmodel::FemModel, solutionstring::Symbol)#{{{
	#independent variable
	α = md.inversion.independent
	#initialize derivative as 0
	∂J_∂α = make_zero(α)
	if md.inversion.onlygrad
		# only compute the gradient
		ComputeGradient!(∂J_∂α, α, femmodel)
		#Put gradient in results
		InputUpdateFromVectorx(femmodel, ∂J_∂α, GradientEnum, VertexSIdEnum)
		RequestedOutputsx(femmodel, [GradientEnum])
	else
		# optimization
		# define cost function and gradient
		# need to build connection between md and x
		f(x) = begin
			fem=DJUICE.ModelProcessor(md, solutionstring)
			DJUICE.CostFunction(x, fem)
		end

		g!(gx, x) = begin
			fem=DJUICE.ModelProcessor(md, solutionstring)
			DJUICE.ComputeGradient!(gx, x, fem)
		end
		nlp = NLPModel(
							α,
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

		independent_enum = StringToEnum(md.inversion.independent_string)
		InputUpdateFromVectorx(femmodel, results_qn.solution, independent_enum, VertexSIdEnum)
		RequestedOutputsx(femmodel, [independent_enum])
	end
end#}}}
function ComputeGradient!(∂J_∂α::Vector{Float64}, α::Vector{Float64}, femmodel::FemModel) #{{{
	# zero ALL depth of the model, make sure we get correct gradient
	dfemmodel = make_zero(Base.Core.Typeof(femmodel), IdDict(), femmodel)
	# zero the gradient
	∂α = make_zero(α)
	# compute the gradient
	autodiff(set_runtime_activity(Enzyme.Reverse), CostFunction, Active, Duplicated(α, ∂α), Duplicated(femmodel,dfemmodel))
	# put gradient back
	∂J_∂α .= ∂α
end#}}}
function CostFunction(α::Vector{Float64}, femmodel::FemModel) #{{{
	# get the md.inversion.independent_string
	control_string = FindParam(String, femmodel.parameters, InversionControlParametersEnum)
	# get the Enum
	controlvar_enum = StringToEnum(control_string)
	if isnothing(controlvar_enum)
		error(control_string, " is not defined in DJUICE, therefore the derivative with respect to ", control_string, " is meaningless")
	end

	# get the cost function list from md.inversion.dependent_string
	cost_list = FindParam(Vector{String}, femmodel.parameters, InversionCostFunctionsEnum)
	cost_enum_list = Vector{IssmEnum}(undef, length(cost_list))
	for (index, value) in enumerate(cost_list)
		cost_enum_list[index] = StringToEnum(value)
	end

	# compute cost function
	# TODO: loop through all controls with respect to all the components in the cost function
	solutionstring = FindParam(Symbol, femmodel.parameters, SolutionTypeEnum)
	# return J
	CostFunctionx(femmodel, α, controlvar_enum, VertexSIdEnum, cost_enum_list, Val(solutionstring))
end#}}}
function CostFunctionx(femmodel::FemModel, α::Vector{Float64}, controlvar_enum::IssmEnum, SId_enum::IssmEnum, cost_enum_list::Vector{IssmEnum}, ::Val{solutionstring}) where solutionstring #{{{
	#Update FemModel accordingly
	InputUpdateFromVectorx(femmodel, α, controlvar_enum, SId_enum)

	#solve PDE
	if(solutionstring===:StressbalanceSolution)
		analysis = StressbalanceAnalysis()
	elseif (solutionstring===:TransientSolution)
		analysis = TransientAnalysis()
	else
		error("not supported")
	end

	Core(analysis, femmodel)

	#update all values of the cost functions
	RequestedOutputsx(femmodel, cost_enum_list)

	#Compute cost function
	J = 0.0
	for i in 1:length(cost_enum_list)
		J += femmodel.results[end-i+1].value
	end

	#return cost function
	return J
end#}}}

