
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

# The misfit functions
function SurfaceAbsVelMisfitx(femmodel::FemModel) #{{{

   #Initialize output
   J = 0.0

   #Sum all element values
   for i in 1:length(femmodel.elements)

      #Get current element
      element = femmodel.elements[i]

      #Should we skip?
      if(!IsIceInElement(femmodel.elements[i])) continue end

      #Retrieve all inputs and parameters
      xyz_list = GetVerticesCoordinates(element.vertices)
      vx_input     = GetInput(element, VxEnum)
      vy_input     = GetInput(element, VyEnum)
      vx_obs_input = GetInput(element, VxObsEnum)
      vy_obs_input = GetInput(element, VyObsEnum)

      #Start integrating
      gauss = GaussTria(3)
      for ig in 1:gauss.numgauss

         Jdet   = JacobianDeterminant(xyz_list, gauss)

         vx    = GetInputValue(vx_input, gauss, ig)
         vy    = GetInputValue(vy_input, gauss, ig)
         vxobs = GetInputValue(vx_obs_input, gauss, ig)
         vyobs = GetInputValue(vy_obs_input, gauss, ig)

         J += gauss.weights[ig]*Jdet*(0.5*(vx-vxobs)^2 + 0.5*(vy-vyobs)^2)
      end
   end

   return J
end#}}}
