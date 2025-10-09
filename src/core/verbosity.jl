global verbositylevel::Int = 0

VerboseMProcessor()::Bool = (verbositylevel & 1)>0
VerboseModule()::Bool = (verbositylevel & 2)>0
VerboseSolution()::Bool = (verbositylevel & 4)>0
VerboseSolver()::Bool = (verbositylevel & 8)>0
VerboseConvergence()::Bool = (verbositylevel & 16)>0
VerboseControl()::Bool = (verbositylevel & 32)>0
VerboseAutodiff()::Bool = (verbositylevel & 128)>0

function SetVerbosityLevel(md::model)
	binary = 0
	if (md.verbose.mprocessor) binary = binary | 1 end
	if (md.verbose.modules) binary = binary | 2 end
	if (md.verbose.solution) binary = binary | 4 end
	if (md.verbose.solver) binary = binary | 8 end
	if (md.verbose.convergence) binary = binary | 16 end
	if (md.verbose.control) binary = binary | 32 end
	if (md.verbose.autodiff) binary = binary | 128 end

	global verbositylevel = binary
	return nothing
end
