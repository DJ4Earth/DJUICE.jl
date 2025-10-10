#Model Processor and Core I/O
function FetchDataToInput(md::model, inputs::Inputs, elements::Vector{Tria}, data::Vector{Float64}, enum::IssmEnum, scaling::Float64=1.) #{{{
	for i in 1:length(elements)
		InputCreate(elements[i], inputs, data, enum, scaling)
	end
	return nothing
end#}}}
function FetchDataToInput(md::model, inputs::Inputs, elements::Vector{Tria}, data::Matrix{Float64}, enum::IssmEnum, scaling::Float64=1.) #{{{
	for i in 1:length(elements)
		InputCreate(elements[i], inputs, data, enum, scaling)
	end
	return nothing
end#}}}
function ModelProcessor(md::model, solutionstring::Symbol) #{{{

	# Set Verbosity once for all
	SetVerbosityLevel(md)

	#Initialize structures
	elements    = Vector{Tria}(undef,0)
	vertices    = Vector{Vertex}(undef,0)
	results     = Vector{Result}(undef,0)
	parameters  = Parameters()
	inputs      = Inputs(md.mesh.numberofelements, md.mesh.numberofvertices, Dict{IssmEnum,Input}())

	#Create  elements, vertices and materials (independent of the analysis)
	CreateElements(elements, md)
	CreateVertices(vertices, md)
	CreateParameters(parameters, md, solutionstring)
	CreateInputs(inputs, elements, md)
	if solutionstring===:TransientSolution
		UpdateParameters(TransientAnalysis(), parameters, md)
	end

	#Now create analysis specific data structure
	if solutionstring===:StressbalanceSolution
		analyses = Analysis[StressbalanceAnalysis()]
	elseif solutionstring===:TransientSolution
		if md.transient.ismovingfront
			analyses = Analysis[StressbalanceAnalysis(), MasstransportAnalysis(), LevelsetAnalysis(), L2ProjectionBaseAnalysis(), ExtrapolationAnalysis()]
		else
			analyses = Analysis[StressbalanceAnalysis(), MasstransportAnalysis()]
		end
	else
		error(solutionstring, " not supported by ModelProcessor")
	end

	#Initialize analysis specific datasets
	numanalyses = length(analyses)
	nodes       = Vector{Vector{Node}}(undef,numanalyses)
	constraints = Vector{Vector{AbstractConstraint}}(undef,numanalyses)
	for i in 1:numanalyses
		analysis = analyses[i]
		if (VerboseMProcessor()) println("   creating datasets for analysis ", typeof(analysis)) end
		nodes[i]       = Vector{Node}(undef,0)
		constraints[i] = Vector{Constraint}(undef,0)

		UpdateParameters(analysis, parameters, md)
		UpdateElements(analysis, elements, inputs, md)
		CreateNodes(analysis, nodes[i], md)
		CreateConstraints(analysis, constraints[i], md)

		#Configure objects
		ConfigureObjectx(elements, nodes[i], vertices, parameters, inputs, i)
	end

	#Inversion?
	if md.inversion.iscontrol
		AddParam(parameters, md.inversion.independent_string, InversionControlParametersEnum)
		AddParam(parameters, md.inversion.dependent_string, InversionCostFunctionsEnum)
		# need to load obs vx and vy
		if issubset(["SurfaceAbsVelMisfit"], md.inversion.dependent_string)
			FetchDataToInput(md, inputs, elements, md.inversion.vx_obs./md.constants.yts,VxObsEnum)
			FetchDataToInput(md, inputs, elements, md.inversion.vy_obs./md.constants.yts,VyObsEnum)
		end
	end

	#Build FemModel
	femmodel = FemModel(analyses, elements, vertices,
							  Vector{Node}(undef,0), nodes,
							  parameters, inputs,
							  Vector{Constraint}(undef,0), constraints,
							  results)

	if (VerboseMProcessor()) println("      detecting active vertices") end
	GetMaskOfIceVerticesLSMx0(femmodel)

	return femmodel

end# }}}
function CreateElements(elements::Vector{Tria},md::model) #{{{

	#Make sure elements is currently empty
	@assert length(elements)==0

	count = 0
	for i in 1:(md.mesh.numberofelements::Int64)

		#Assume Linear Elements for now
		vertexids = (md.mesh.elements[i,:]::Vector{Int64})
		nodeids   = (md.mesh.elements[i,:]::Vector{Int64})

		#Call constructor and add to dataset elements
		push!(elements,Tria(i,count, vertexids))
	end

	return nothing
end# }}}
function CreateVertices(vertices::Vector{Vertex},md::model) #{{{

	#Make sure vertices is currently empty
	@assert length(vertices)==0

	#Get data from md
	x = md.mesh.x
	y = md.mesh.y

	count = 0
	for i in 1:md.mesh.numberofvertices
		push!(vertices,Vertex(i,x[i],y[i],0.))
	end

	return nothing
end# }}}
function CreateParameters(parameters::Parameters,md::model,solutionstring::Symbol) #{{{

	#Get data from md
	AddParam(parameters,md.materials.rho_ice,MaterialsRhoIceEnum)
	AddParam(parameters,md.materials.rho_water,MaterialsRhoSeawaterEnum)
	AddParam(parameters,md.materials.rho_freshwater,MaterialsRhoFreshwaterEnum)
	AddParam(parameters,md.constants.g,ConstantsGEnum)

	#some parameters that did not come with the iomodel
	AddParam(parameters, solutionstring, SolutionTypeEnum)
	
	#Set step and time, this will be overwritten if we run a transient
	AddParam(parameters,1,StepEnum)
	AddParam(parameters,0.0,TimeEnum)

	#Is moving front 
	AddParam(parameters,md.transient.ismovingfront,TransientIsmovingfrontEnum)
	#Is grounding line
	AddParam(parameters,md.transient.isgroundingline,TransientIsgroundinglineEnum)
	AddParam(parameters,md.groundingline.migration,GroundinglineMigrationEnum)

	return nothing
end# }}}
function CreateInputs(inputs::Inputs,elements::Vector{Tria},md::model) #{{{

	#Only assume we have Matice for now
	FetchDataToInput(md,inputs,elements,md.materials.rheology_B,MaterialsRheologyBEnum)
	FetchDataToInput(md,inputs,elements,md.materials.rheology_n,MaterialsRheologyNEnum)

	return nothing
end# }}}
function OutputResultsx(femmodel::FemModel, md::model, solutionkey::Symbol)# {{{

	if solutionkey===:TransientSolution

		#Compute maximum number of steps
		maxstep = 0
		for i in length(femmodel.results)
			if(femmodel.results[i].step>maxstep) maxstep = femmodel.results[i].step end
		end

		#Initialize vector now that we know the size
		output = Vector{Dict}(undef, maxstep)
		for i in 1:maxstep; output[i] = Dict() end

		#Insert results in vector
		for i in 1:length(femmodel.results)
			result = femmodel.results[i]
			step   = femmodel.results[i].step
			scale = OutputEnumToScale(md, result.enum)
			(output[step])[EnumToString(result.enum)] = result.value .* scale
		end
	else
		output = Dict()
		for i in 1:length(femmodel.results)
			result = femmodel.results[i]
			scale = OutputEnumToScale(md, result.enum)
			output[EnumToString(result.enum)] = result.value .* scale
		end
	end

	md.results[string(solutionkey)] = output

	return nothing
end# }}}

#Other modules
function ConfigureObjectx(elements::Vector{Tria}, nodes::Vector{Node}, vertices::Vector{Vertex}, parameters::Parameters, inputs::Inputs, analysis::Int64) #{{{

	#Configure elements
	for i in 1:length(elements)
		Configure(elements[i], nodes, vertices, parameters, inputs, analysis)
	end

	#Configure inputs
	Configure(inputs, parameters)

	return nothing
end# }}}
function CreateNodalConstraintsx(nodes::Vector{Node})# {{{

	#Allocate vector
	ssize=NumberOfDofs(nodes,SsetEnum)
	ys=IssmVector(ssize)

	#constraints vector with the constraint values
	for i in 1:length(nodes)
		CreateNodalConstraints(nodes[i],ys)
	end

	return ys
end# }}}
function FloatingiceMeltingRate(femmodel::FemModel)# {{{

	basalforcing_model = FindParam(IssmEnum, femmodel.parameters, BasalforcingsEnum)

	if basalforcing_model==FloatingMeltRateEnum
		#Nothing to be done
	elseif basalforcing_model==LinearFloatingMeltRateEnum
		LinearFloatingiceMeltingRate(femmodel)
	end

end# }}}
function LinearFloatingiceMeltingRate(femmodel::FemModel)# {{{

	values       = Vector{Float64}(undef,3)
	base         = Vector{Float64}(undef,3)
	perturbation = Vector{Float64}(undef,3)

	for i=1:length(femmodel.elements)
		element = femmodel.elements[i]

		deepwaterel    = FindParam(Float64, element, BasalforcingsDeepwaterElevationEnum)
		deepwatermelt  = FindParam(Float64, element, BasalforcingsDeepwaterMeltingRateEnum)
		upperwaterel   = FindParam(Float64, element, BasalforcingsUpperwaterElevationEnum)
		upperwatermelt = FindParam(Float64, element, BasalforcingsUpperwaterMeltingRateEnum)

		GetInputListOnVertices!(element, base, BaseEnum)
		GetInputListOnVertices!(element, perturbation, BasalforcingsPerturbationMeltingRateEnum)
		for i in 1:3
			if(base[i]>=upperwaterel)
				values[i]=upperwatermelt
			elseif (base[i]<deepwaterel)
				values[i]=deepwatermelt
			else
				alpha = (base[i]-upperwaterel)/(deepwaterel-upperwaterel)
				values[i]=deepwatermelt*alpha+(1.0 - alpha)*upperwatermelt
			end

			#Add perturbation
			values[i] += perturbation[i]
		end

		AddInput(element, BasalforcingsFloatingiceMeltingRateEnum, values, P1Enum)
	end

end# }}}
function GetMaskOfIceVerticesLSMx0(femmodel::FemModel) #{{{

	#Initialize vector with number of vertices
	numvertices=length(femmodel.vertices)
	if(numvertices==0) return end

	#Initialize vector
	nbv = 3
	onesvec      = ones(nbv)
	vec_mask_ice = IssmVector(numvertices)
	vertexids    = Vector{Int64}(undef, nbv)

	#Assign values to vector
	for i in 1:length(femmodel.elements)
		if (IsIceInElement(femmodel.elements[i]))
			for j in 1:nbv
				vertexids[j] = femmodel.elements[i].vertices[j].sid
			end
			SetValues!(vec_mask_ice, nbv, vertexids, onesvec)
		end
	end
	Assemble!(vec_mask_ice)

	#Serialize vector
	vec_mask_ice_serial = ToSerial(vec_mask_ice)

	#Update IceMaskNodeActivationEnum in elements
	InputUpdateFromVectorx(femmodel, vec_mask_ice_serial, IceMaskNodeActivationEnum, VertexSIdEnum)

	return nothing
end#}}}
function GetSolutionFromInputsx(analysis::Analysis,femmodel::FemModel) #{{{

	#Get size of vector
	gsize = NumberOfDofs(femmodel.nodes,GsetEnum)

	#Initialize solution vector
	ug = IssmVector(gsize)

	#Go through elements and plug in solution
	for i=1:length(femmodel.elements)
		GetSolutionFromInputs(analysis,ug,femmodel.elements[i])
	end

	return ug
end#}}}
function IceVolumeAboveFloatationx(femmodel::FemModel)# {{{

	total_ice_volume_af = 0.0
	for i=1:length(femmodel.elements)
		total_ice_volume_af += IceVolumeAboveFloatation(femmodel.elements[i])
	end
	return total_ice_volume_af
end# }}}
function IceVolumex(femmodel::FemModel)# {{{

	total_ice_volume = 0.0
	for i=1:length(femmodel.elements)
		total_ice_volume += IceVolume(femmodel.elements[i])
	end
	return total_ice_volume
end# }}}
function InputUpdateFromSolutionx(analysis::Analysis,ug::IssmVector,femmodel::FemModel) #{{{

	#Go through elements and plug in solution
	for i=1:length(femmodel.elements)
		InputUpdateFromSolution(analysis,ug.vector,femmodel.elements[i])
	end

	return Nothing
end#}}}
function InputUpdateFromVectorx(femmodel::FemModel, vector::Vector{Float64}, enum::IssmEnum, layout::IssmEnum)# {{{

	#Go through elements and plug in solution
	for i=1:length(femmodel.elements)
		InputUpdateFromVector(femmodel.elements[i], vector, enum, layout)
	end

	return nothing
end#}}}
function InputDuplicatex(femmodel::FemModel, oldenum::IssmEnum, newenum::IssmEnum) #{{{
	DuplicateInput(femmodel.inputs, oldenum, newenum)
	return nothing
end#}}}
function Mergesolutionfromftogx(ug::IssmVector, uf::IssmVector, ys::IssmVector, nodes::Vector{Node}) #{{{

	#Go through elements and plug in solution
	for i=1:length(nodes)
		VecMerge(nodes[i],ug,uf.vector,ys.vector)
	end

	return ug

end#}}}
function MigrateGroundinglinex(femmodel::FemModel)# {{{

	for i=1:length(femmodel.elements)
		MigrateGroundingLine(femmodel.elements[i])
	end

	return nothing
end# }}}
function NodesDofx(nodes::Vector{Node}, parameters::Parameters) #{{{

	#Do we have any nodes?
	if length(nodes)==0
		return
	end

	#Do we really need to update dof indexing
	if(~RequiresDofReindexing(nodes)) return end

	if (VerboseModule()) print("   Renumbering degrees of freedom\n") end
	DistributeDofs(nodes,GsetEnum)
	DistributeDofs(nodes,FsetEnum)
	DistributeDofs(nodes,SsetEnum)

	return nothing
end# }}}
function OutputEnumToScale(md::model, result::IssmEnum) #{{{
	if (result == VxEnum) scale = md.constants.yts
	elseif (result == VyEnum) scale = md.constants.yts
	elseif (result == VzEnum) scale = md.constants.yts
	elseif (result == VelEnum) scale = md.constants.yts
	elseif (result == CalvingCalvingrateEnum) scale = md.constants.yts
	else scale = 1.0
	end

	return scale
end#}}}
function Reduceloadx!(pf::IssmVector, Kfs::IssmMatrix, ys::IssmVector) #{{{

	#Is there anything to do?
	m, n = GetSize(Kfs)

	if(m*n>0)

		#Allocate Kfs*ys
		Kfsy_s=IssmVector(m)

		#Perform multiplication
		MatMult!(Kfs,ys,Kfsy_s)

		#Subtract Kfs*ys from pf
		AXPY!(pf,-1.0,Kfsy_s)

	end
	return nothing
end#}}}
function Reducevectorgtofx(ug::IssmVector,nodes::Vector{Node}) #{{{

	#Get size of output vector
	fsize = NumberOfDofs(nodes,FsetEnum)

	#Initialize output vector
	uf = IssmVector(fsize)

	#Go through elements and plug in solution
	for i=1:length(nodes)
		VecReduce(nodes[i],ug.vector,uf)
	end

	return uf

end#}}}
function RequestedOutputsx(femmodel::FemModel,outputlist::Vector{IssmEnum})# {{{

	#Get Step and Time from parameters
	step = FindParam(Int64, femmodel.parameters,StepEnum)
	time = FindParam(Float64, femmodel.parameters,TimeEnum)


	#Now fetch results
	for i in 1:length(outputlist)
		#result is a vector
		isvec = false
		#See if outputlist[i] is an input
		if outputlist[i] == IceVolumeEnum
			double_result = IceVolumex(femmodel)
		elseif outputlist[i] == IceVolumeAboveFloatationEnum
			double_result = IceVolumeAboveFloatationx(femmodel)
		elseif outputlist[i] == SurfaceAbsVelMisfitEnum
			double_result = SurfaceAbsVelMisfitx(femmodel)
		elseif outputlist[i] == DragCoefficientAbsGradientEnum #TODO: define a new Enum
			double_result = ControlVariableAbsGradientx(femmodel)
		elseif outputlist[i]>InputsSTARTEnum && outputlist[i]<InputsENDEnum # default
			#Create Result
			input  = GetInput(femmodel.inputs, outputlist[i])
			input_copy = Vector{Float64}(undef, length(input.values))
			copyto!(input_copy, input.values)
			result = Result(outputlist[i], step, time, input_copy)
			# flag of vector
			isvec = true
		else
			error("Output ",outputlist[i]," not supported yet")
		end

		if (!isvec)
			result = Result(outputlist[i], step, time, double_result)
		end
		#Add to femmodel.results dataset
		AddResult!(femmodel, result)
		#		push!(femmodel.results, result)
	end
	return nothing
end# }}}
function SetActiveNodesLSMx(femmodel::FemModel) #{{{

	#Check mask of each element to see if element is active
	numnodes = 3
	mask = Vector{Float64}(undef, numnodes)
	for i in 1:length(femmodel.elements)
		GetInputListOnNodes!(femmodel.elements[i], mask, IceMaskNodeActivationEnum)
		for j in 1:numnodes
			node = femmodel.elements[i].nodes[j]
			if(mask[j]==1.) Activate!(node)
			else             Deactivate!(node)
			end
		end
	end
	return nothing
end#}}}
function SpcNodesx(nodes::Vector{Node},constraints::Vector{AbstractConstraint},parameters::Parameters) #{{{

	for i in 1:length(constraints)
		ConstrainNode(constraints[i],nodes,parameters)
	end

	return nothing
end# }}}
function SystemMatricesx(femmodel::FemModel,analysis::Analysis)# {{{

	#Allocate matrices
	if (VerboseModule()) println("   Allocating matrices") end
	fsize = NumberOfDofs(femmodel.nodes,FsetEnum)
	ssize = NumberOfDofs(femmodel.nodes,SsetEnum)
	Kff = IssmMatrix(fsize,fsize)
	Kfs = IssmMatrix(fsize,ssize)
	pf  = IssmVector(fsize)

	#Construct Stiffness matrix and load vector from elements
	if (VerboseModule()) println("   Assembling matrices") end
	# Determine the type of analysis, decide if override icemask 
	if typeof(analysis) == typeof(ExtrapolationAnalysis())
		overrideicemask = true
	elseif typeof(analysis) == typeof(LevelsetAnalysis())
		overrideicemask = true
	elseif typeof(analysis) == typeof(L2ProjectionBaseAnalysis())
		overrideicemask = true
	else
		overrideicemask = false
	end


	for i in 1:length(femmodel.elements)

		if (IsIceInElement(femmodel.elements[i]) || overrideicemask)
			if (AnyFSet(femmodel.elements[i]) || (typeof(analysis) == typeof(StressbalanceAnalysis())))
				Ke = CreateKMatrix(analysis,femmodel.elements[i])
				if Ke != nothing
					AddToGlobal!(Ke,Kff,Kfs)
				end

				pe = CreatePVector(analysis,femmodel.elements[i])
				if pe != nothing
					AddToGlobal!(pe,pf)
				end
			end
		end
	end

	Assemble!(Kff)
	Assemble!(Kfs)
	Assemble!(pf)

	return Kff, Kfs, pf
end# }}}
function UpdateConstraintsx(femmodel::FemModel, analysis::Analysis)# {{{

	#First, see if the analysis needs to change constraints
	UpdateConstraints(analysis, femmodel)

	#Second, constraints might be time dependent
	SpcNodesx(femmodel.nodes, femmodel.constraints, femmodel.parameters)

	#Now, update degrees of freedoms
	NodesDofx(femmodel.nodes, femmodel.parameters)

	return nothing
end# }}}
