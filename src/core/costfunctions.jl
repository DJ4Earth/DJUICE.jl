
function costfunction(femmodel::FemModel, α::Vector{Float64}) #{{{
	# TODO: automatic determine costfunction from md.inversion
	controlstring = FindParam(String, femmodel.parameters, InversionControlParametersEnum)
	if (controlstring == "FrictionC")
		controlEnum = FrictionCoefficientEnum
	elseif (controlstring == "RheologyB")
		controlEnum = MaterialsRheologyBEnum
	else
		error(controlstring, " is not supported, just for now. ")
	end
	costfunction(femmodel, α,  controlEnum, VertexSIdEnum)
end#}}}
function costfunction(femmodel::FemModel, α::Vector{Float64}, variableEnum::IssmEnum, SIdEnum::IssmEnum) #{{{
   #Update FemModel accordingly
	InputUpdateFromVectorx(femmodel, α, variableEnum, SIdEnum)

   #solve PDE
   analysis = StressbalanceAnalysis()
   Core(analysis, femmodel)

   #Compute cost function
   J = SurfaceAbsVelMisfitx(femmodel)

   #return cost function
   return J
end#}}}
