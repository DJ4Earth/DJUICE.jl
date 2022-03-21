#Model Processor and Core I/O
function FetchDataToInput(md::model,inputs::Inputs,elements::Vector{Tria},data::Vector{Float64},enum::IssmEnum) #{{{
	for i in 1:length(elements)
		InputCreate(elements[i],inputs,data,enum)
	end
	return nothing
end#}}}
function ModelProcessor(md::model, solutionstring::String) #{{{

	#Initialize structures
	elements    = Vector{Tria}(undef,0)
	vertices    = Vector{Vertex}(undef,0)
	results     = Vector{Result}(undef,0)
	parameters  = Parameters()
	inputs      = Inputs(md.mesh.numberofelements, md.mesh.numberofvertices, Dict{IssmEnum,Input}())

	#Create  elements, vertices and materials (independent of the analysis)
	CreateElements(elements, md)
	CreateVertices(vertices, md)
	CreateParameters(parameters, md)
	CreateInputs(inputs,elements, md)
	if solutionstring=="TransientSolution"
		UpdateParameters(TransientAnalysis(), parameters, md)
	end

	#Now create analysis specific data structure
	if solutionstring=="StressbalanceSolution"
		analyses = [StressbalanceAnalysis()]
	elseif solutionstring=="TransientSolution"
		analyses = [StressbalanceAnalysis(), MasstransportAnalysis()]
	else
		error(solutionstring, " not supported by ModelProcessor")
	end

	#Initialize analysis specific datasets
	numanalyses = length(analyses)
	nodes       = Vector{Vector{Node}}(undef,numanalyses)
	constraints = Vector{Vector{Constraint}}(undef,numanalyses)
	for i in 1:numanalyses
		analysis = analyses[i]
		println("   creating datasets for analysis ", typeof(analysis))
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
		FetchDataToInput(md, inputs, elements, md.inversion.vx_obs./md.constants.yts,VxObsEnum)
		FetchDataToInput(md, inputs, elements, md.inversion.vy_obs./md.constants.yts,VyObsEnum)
	end

	#Build FemModel
	femmodel = FemModel(analyses, elements, vertices,
							  Vector{Node}(undef,0), nodes,
							  parameters, inputs,
							  Vector{Constraint}(undef,0), constraints,
							  results)

	println("      detecting active vertices")
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
function CreateParameters(parameters::Parameters,md::model) #{{{

	#Get data from md
	AddParam(parameters,md.materials.rho_ice,MaterialsRhoIceEnum)
	AddParam(parameters,md.materials.rho_water,MaterialsRhoSeawaterEnum)
	AddParam(parameters,md.materials.rho_freshwater,MaterialsRhoFreshwaterEnum)
	AddParam(parameters,md.constants.g,ConstantsGEnum)
	
	#Set step and time, this will be overwritten if we run a transient
	AddParam(parameters,1,StepEnum)
	AddParam(parameters,0.0,TimeEnum)

	return nothing
end# }}}
function CreateInputs(inputs::Inputs,elements::Vector{Tria},md::model) #{{{

	#Only assume we have Matice for now
	FetchDataToInput(md,inputs,elements,md.materials.rheology_B,MaterialsRheologyBEnum)
	FetchDataToInput(md,inputs,elements,md.materials.rheology_n,MaterialsRheologyNEnum)

	return nothing
end# }}}
function OutputResultsx(femmodel::FemModel, md::model, solution::String)# {{{


	if solution=="TransientSolution"

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
			(output[step])[EnumToString(result.enum)] = result.value
		end
	else
		output = Dict()
		for i in length(femmodel.results)
			result = femmodel.results[i]
			output[EnumToString(result.enum)] = result.value
		end
	end

	md.results[solution] = output

	return nothing
end# }}}

#Other modules
function ConfigureObjectx(elements::Vector{Tria}, nodes::Vector{Node}, vertices::Vector{Vertex}, parameters::Parameters, inputs::Inputs, analysis::Int64) #{{{

	for i in 1:length(elements)
		Configure(elements[i], nodes, vertices, parameters, inputs, analysis)
	end

	return nothing
end# }}}
function SpcNodesx(nodes::Vector{Node},constraints::Vector{Constraint},parameters::Parameters) #{{{

	for i in 1:length(constraints)
		ConstrainNode(constraints[i],nodes,parameters)
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

	print("   Renumbering degrees of freedom\n")
	DistributeDofs(nodes,GsetEnum)
	DistributeDofs(nodes,FsetEnum)
	DistributeDofs(nodes,SsetEnum)

	return nothing
end# }}}
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
function InputUpdateFromSolutionx(analysis::Analysis,ug::IssmVector,femmodel::FemModel) #{{{

	#Go through elements and plug in solution
	for i=1:length(femmodel.elements)
		InputUpdateFromSolution(analysis,ug.vector,femmodel.elements[i])
	end

	return ug

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
function Mergesolutionfromftogx(ug::IssmVector, uf::IssmVector, ys::IssmVector, nodes::Vector{Node}) #{{{

	#Go through elements and plug in solution
	for i=1:length(nodes)
		VecMerge(nodes[i],ug,uf.vector,ys.vector)
	end

	return ug

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
function SystemMatricesx(femmodel::FemModel,analysis::Analysis)# {{{

	#Allocate matrices
	fsize = NumberOfDofs(femmodel.nodes,FsetEnum)
	ssize = NumberOfDofs(femmodel.nodes,SsetEnum)
	Kff = IssmMatrix(fsize,fsize)
	Kfs = IssmMatrix(fsize,ssize)
	pf  = IssmVector(fsize)

	#Construct Stiffness matrix and load vector from elements
	for i in 1:length(femmodel.elements)
		Ke = CreateKMatrix(analysis,femmodel.elements[i])
		pe = CreatePVector(analysis,femmodel.elements[i])

		if(!isnothing(Ke)) AddToGlobal!(Ke,Kff,Kfs) end
		if(!isnothing(pe)) AddToGlobal!(pe,pf) end
	end

	Assemble!(Kff)
	Assemble!(Kfs)
	Assemble!(pf)
	
	return Kff, Kfs, pf
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
function RequestedOutputsx(femmodel::FemModel,outputlist::Vector{IssmEnum})# {{{

	#Get Step and Time from parameters
	step = FindParam(Int64, femmodel.parameters,StepEnum)
	time = FindParam(Float64, femmodel.parameters,TimeEnum)
	
	#Now fetch results
	for i in 1:length(outputlist)

		#See if outputlist[i] is an input
		if outputlist[i]>InputsSTARTEnum && outputlist[i]<InputsENDEnum

			#Create Result
			input  = GetInput(femmodel.inputs,outputlist[i])
			result = Result(outputlist[i], step, time, copy(input.values))

			#Add to femmodel.results dataset
			push!(femmodel.results, result)

		else
			println("Output ",outputlist[i]," not supported yet")
		end
	end
	return nothing
end# }}}
function MigrateGroundinglinex(femmodel::FemModel)# {{{

	for i=1:length(femmodel.elements)
		MigrateGroundingLine(femmodel.elements[i])
	end

	return nothing
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
