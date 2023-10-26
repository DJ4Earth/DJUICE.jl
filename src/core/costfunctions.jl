
function costfunction(femmodel::FemModel, α::Vector{Float64}) #{{{
	# TODO: automatic determine costfunction from md.inversion
	costfunction(femmodel, α,  FrictionCoefficientEnum, VertexSIdEnum)
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
