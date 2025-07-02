#L2ProjectionBaseAnalysis class definition
struct L2ProjectionBaseAnalysis <: Analysis#{{{
end #}}}

#Model Processing
function CreateConstraints(analysis::L2ProjectionBaseAnalysis,constraints::Vector{AbstractConstraint},md::model) #{{{
	return nothing
end#}}}
function CreateNodes(analysis::L2ProjectionBaseAnalysis,nodes::Vector{Node},md::model) #{{{

	numdof = 1
	for i in 1:md.mesh.numberofvertices
		push!(nodes,Node(i,i,true,true,numdof,-ones(Int64,numdof), ones(Int64,numdof), -ones(Int64,numdof), zeros(numdof)))
	end

	return nothing
end#}}}
function UpdateElements(analysis::L2ProjectionBaseAnalysis,elements::Vector{Tria}, inputs::Inputs, md::model) #{{{

	#Provide node indices to element
	for i in 1:md.mesh.numberofelements
		Update(elements[i],inputs,i,md,P1Enum)
	end

	#Add necessary inputs to perform this analysis
	FetchDataToInput(md,inputs,elements,md.geometry.surface,SurfaceEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.base,BaseEnum)
	FetchDataToInput(md,inputs,elements,md.mask.ice_levelset, MaskIceLevelsetEnum)

	return nothing
end#}}}
function UpdateParameters(analysis::L2ProjectionBaseAnalysis,parameters::Parameters,md::model) #{{{
	return nothing
end#}}}

#Finite Element Analysis
function Core(analysis::L2ProjectionBaseAnalysis,femmodel::FemModel)# {{{
	return nothing
end #}}}
function CreateKMatrix(analysis::L2ProjectionBaseAnalysis,element::Tria)# {{{
	#Internmediaries
	numnodes = 3
	
	#Initialize Element matrix and basis function derivatives
	Ke = ElementMatrix(element.nodes)
	basis  = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)

	#Start integrating
	gauss = GaussTria(2)
	for ig in 1:gauss.numgauss

		Jdet = JacobianDeterminant(xyz_list, gauss)
		NodalFunctions(element, basis, gauss, ig, P1Enum)

      #Transient term
		for i in 1:numnodes
			for j in 1:numnodes
				Ke.values[i ,j] += gauss.weights[ig]*Jdet*basis[i]*basis[j]
			end
		end
	end

	return Ke
end #}}}
function CreatePVector(analysis::L2ProjectionBaseAnalysis,element::Tria)# {{{

	#Internmediaries
	numnodes = 3

	#Initialize Element vectro and basis functions
	pe = ElementVector(element.nodes)
	basis = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	input_enum = FindParam(IssmEnum, element, InputToL2ProjectEnum)

	if input_enum == SurfaceSlopeXEnum || input_enum == SurfaceSlopeYEnum
		input2 = GetInput(element, SurfaceEnum)
	elseif input_enum == BedSlopeXEnum || input_enum == BedSlopeYEnum || input_enum == BaseSlopeXEnum || input_enum == BaseSlopeYEnum
		input2 = GetInput(element, BaseEnum)
	elseif input_enum == LevelsetfunctionSlopeXEnum || input_enum == LevelsetfunctionSlopeYEnum
		input2 = GetInput(element, MaskIceLevelsetEnum)
	else
		input2 = GetInput(element, input_enum)
	end

	#Start integrating
	gauss = GaussTria(2)
	for ig in 1:gauss.numgauss
		Jdet = JacobianDeterminant(xyz_list, gauss)
		#Get nodal basis
		NodalFunctions(element, basis, gauss, ig, P1Enum)
		
      slopes = GetInputDerivativeValue(input2, xyz_list, gauss, ig)
		if input_enum == SurfaceSlopeXEnum || input_enum == BedSlopeXEnum || input_enum == BaseSlopeXEnum || input_enum == LevelsetfunctionSlopeXEnum
			value = slopes[1]
		elseif input_enum == SurfaceSlopeXEnum || input_enum == BedSlopeXEnum || input_enum == BaseSlopeXEnum || input_enum == LevelsetfunctionSlopeXEnum
			value = slopes[2]
		else
			value = GetInputValue(input2, gauss, ig)
		end

		for i in 1:numnodes
         pe.values[i] += gauss.weights[ig]*Jdet*value*basis[i]
		end
	end

	return pe
end #}}}
function GetSolutionFromInputs(analysis::L2ProjectionBaseAnalysis,ug::IssmVector,element::Tria) #{{{

	#Get dofs for this finite element
	doflist = GetDofList(element,GsetEnum)
	@assert length(doflist)==3

	#Fetch inputs
	mask_input = GetInput(element, MaskIceL2ProjectionBaseEnum)

	#Loop over each node and enter solution in ug
	count = 0
	gauss=GaussTria(P1Enum)
	for i in 1:gauss.numgauss
		mask = GetInputValue(mask_input, gauss, i)

		count += 1
		ug.vector[doflist[count]] = mask
	end

	#Make sure we reached all the values
	@assert count==length(doflist)

	return nothing
end#}}}
function InputUpdateFromSolution(analysis::L2ProjectionBaseAnalysis,ug::Vector{Float64},element::Tria) #{{{
	input_enum = FindParam(IssmEnum, element, InputToL2ProjectEnum)
	InputUpdateFromSolutionOneDof(element, ug, input_enum)
end#}}}
function UpdateConstraints(analysis::L2ProjectionBaseAnalysis, femmodel::FemModel) #{{{
	#Default do nothing
	return nothing
end#}}}
