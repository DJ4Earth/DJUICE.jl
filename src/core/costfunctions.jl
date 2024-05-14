
function costfunction(α::Vector{Float64}, femmodel::FemModel) #{{{
	# get the md.inversion.control_string
	control_string = FindParam(String, femmodel.parameters, InversionControlParametersEnum)
	# get the Enum
	controlvar_enum = StringToEnum(control_string)
	if isnothing(controlvar_enum)
		error(control_string, " is not defined in DJUICE, therefore the derivative with respect to ", control_string, " is meaningless")
	end
	# compute cost function
	# TODO: loop through all controls with respect to all the components in the cost function
	costfunction(α, femmodel, controlvar_enum, VertexSIdEnum)
end#}}}
function costfunction(α::Vector{Float64}, femmodel::FemModel, controlvar_enum::IssmEnum, SId_enum::IssmEnum) #{{{
	#Update FemModel accordingly
	InputUpdateFromVectorx(femmodel, α, controlvar_enum, SId_enum)

	#solve PDE
	analysis = StressbalanceAnalysis()
	Core(analysis, femmodel)

	#Compute cost function
	J = SurfaceAbsVelMisfitx(femmodel)

	J += ControlVariableAbsGradientx(femmodel, α, controlvar_enum)

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
function ControlVariableAbsGradientx(femmodel::FemModel, α::Vector{Float64}, controlvar_enum::IssmEnum) #{{{

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
		controlvar_input = GetInput(element, controlvar_enum)

		#Start integrating
		gauss = GaussTria(3)
		for ig in 1:gauss.numgauss

			Jdet   = JacobianDeterminant(xyz_list, gauss)
			# TODO: add weights
			dα = GetInputDerivativeValue(controlvar_input, xyz_list, gauss, ig)
			J += gauss.weights[ig]*0.5*Jdet*sum(dα.*dα)
		end
	end

	return J
end#}}}
