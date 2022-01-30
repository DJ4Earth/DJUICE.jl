include("../usr/classes.jl")
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

function IssmCore(md::model, solutionstring::String) #{{{

	#Construct FemModel
	femmodel=ModelProcessor(md, solutionstring)

	#Solve FIXME: to be improved later...
	if(solutionstring=="StressbalanceSolution")
		analysis = StressbalanceAnalysis()
	elseif (solutionstring=="TransientSolution")
		analysis = TransientAnalysis()
	else
		error("not supported")
	end
	Core(analysis, femmodel)

	#move results to md
	OutputResultsx(femmodel, md, solutionstring)

end# }}}
