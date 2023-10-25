using Enzyme
Enzyme.API.looseTypeAnalysis!(false)
Enzyme.API.strictAliasing!(false)
Enzyme.API.typeWarning!(false)

include("./issmenums.jl")
include("./toolkits.jl")
include("./gauss.jl")
include("./parameters.jl")
include("./inputs.jl")
include("./vertices.jl")
include("./nodes.jl")
include("./elements.jl")
include("./constraints.jl")
include("./results.jl")
include("./matice.jl")
include("./friction.jl")
include("./analyses/analysis.jl")
include("./femmodel.jl")

#All analyses
include("./analyses/stressbalanceanalysis.jl")
include("./analyses/masstransportanalysis.jl")
include("./analyses/transientanalysis.jl")

include("./solutionsequences.jl")
include("./modules.jl")
include("./elementmatrix.jl")
include("./utils.jl")

function solve(md::model, solution::Symbol) #{{{

	#Process incoming string
	if solution===:sb || solution===:Stressbalance
		solutionkey = :StressbalanceSolution
	elseif solution===:tr || solution===:Transient
		solutionkey = :TransientSolution
	else
		error("solutionkey "*solution*" not supported!");
	end

	#Construct FemModel from md
	femmodel=ModelProcessor(md, solutionkey)

	#Solve (FIXME: to be improved later...)
	if(solutionkey===:StressbalanceSolution)
		analysis = StressbalanceAnalysis()
	elseif (solutionkey===:TransientSolution)
		analysis = TransientAnalysis()
	else
		error("not supported")
	end
	Core(analysis, femmodel)

	#add results to md.results
	OutputResultsx(femmodel, md, solutionkey)

	return md
end# }}}

#Automatic differentiation
function costfunction(femmodel::FemModel, α::Vector{Float64})

	#Update FemModel accordingly
	InputUpdateFromVectorx(femmodel, α, FrictionCoefficientEnum, VertexSIdEnum)

	#solve PDE
	analysis = StressbalanceAnalysis()
	Core(analysis, femmodel)

	#Compute cost function
	J = SurfaceAbsVelMisfitx(femmodel)

	#return cost function
	return J
end
function solve2(md::model,isAD::Bool) #{{{

	#Construct FemModel from md
	femmodel=ModelProcessor(md, :StressbalanceSolution)

	if(isAD)
		#Active variable
		α = md.friction.coefficient
		#initialize derivative as 0
		∂J_∂α = zero(α)

		#Misc Enzyme options
		println("CALLING AUTODIFF, prepare to die...")
		dfemmodel = deepcopy(femmodel)
		@time autodiff(Enzyme.Reverse, costfunction, Duplicated(femmodel, dfemmodel), Duplicated(α, ∂J_∂α))

		#Put gradient in results
		InputUpdateFromVectorx(femmodel, ∂J_∂α, GradientEnum, VertexSIdEnum)
		RequestedOutputsx(femmodel, [GradientEnum])

	else
		analysis = StressbalanceAnalysis()
		Core(analysis, femmodel)
	end

	OutputResultsx(femmodel, md, :StressbalanceSolution)

	return md
end# }}}
