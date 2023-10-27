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
include("./costfunctions.jl")

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
	if (md.inversion.iscontrol) # solve inverse problem
		#independent variable
		controlstring = md.inversion.independent
		if (controlstring == "FrictionC")
			α = md.friction.coefficient
		elseif (controlstring == "RheologyB")
			α = md.materials.rheology_B
		else
			error(controlstring, " is not supported, just for now. ")
		end
		#initialize derivative as 0
		∂J_∂α = zero(α)

		# zero ALL depth of the model, make sure we get correct gradient
		dfemmodel = Enzyme.Compiler.make_zero(Base.Core.Typeof(femmodel), IdDict(), femmodel)
		# compute the gradient
		println("CALLING AUTODIFF, prepare to die...")
		@time autodiff(Enzyme.Reverse, costfunction, Duplicated(femmodel, dfemmodel), Duplicated(α, ∂J_∂α))

		#Put gradient in results
		InputUpdateFromVectorx(femmodel, ∂J_∂α, GradientEnum, VertexSIdEnum)
		RequestedOutputsx(femmodel, [GradientEnum])
	else # otherwise forward problem
		if(solutionkey===:StressbalanceSolution)
			analysis = StressbalanceAnalysis()
		elseif (solutionkey===:TransientSolution)
			analysis = TransientAnalysis()
		else
			error("not supported")
		end
		Core(analysis, femmodel)
	end

	#add results to md.results
	OutputResultsx(femmodel, md, solutionkey)

	return md
end# }}}
