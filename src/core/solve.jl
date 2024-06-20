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
include("./control.jl")

#All analyses
include("./analyses/stressbalanceanalysis.jl")
include("./analyses/masstransportanalysis.jl")
include("./analyses/transientanalysis.jl")
include("./analyses/levelsetanalysis.jl")

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
		Control_Core(md, femmodel)
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
