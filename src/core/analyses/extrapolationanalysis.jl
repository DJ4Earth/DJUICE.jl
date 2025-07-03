#ExtrapolationAnalysis class definition
struct ExtrapolationAnalysis <: Analysis#{{{
end #}}}

#Model Processing
function CreateConstraints(analysis::ExtrapolationAnalysis,constraints::Vector{AbstractConstraint},md::model) #{{{
	# doing nothing for now
	return nothing
end#}}}
function CreateNodes(analysis::ExtrapolationAnalysis,nodes::Vector{Node},md::model) #{{{

	numdof = 1
	for i in 1:md.mesh.numberofvertices
		push!(nodes,Node(i,i,true,true,numdof,-ones(Int64,numdof), ones(Int64,numdof), -ones(Int64,numdof), zeros(numdof)))
	end

	return nothing
end#}}}
function UpdateElements(analysis::ExtrapolationAnalysis,elements::Vector{Tria}, inputs::Inputs, md::model) #{{{

	#Provide node indices to element
	for i in 1:md.mesh.numberofelements
		Update(elements[i],inputs,i,md,P1Enum)
	end

	return nothing
end#}}}
function UpdateParameters(analysis::ExtrapolationAnalysis,parameters::Parameters,md::model) #{{{
	return nothing
end#}}}

#Finite Element Analysis
function Core(analysis::ExtrapolationAnalysis,femmodel::FemModel)# {{{

   SetCurrentConfiguration!(femmodel, analysis)

   #Fetch parameters relevant to solution sequence
   extvar_num  = FindParam(IssmEnum, femmodel.parameters,ExtrapolationVariableEnum)

	SetCurrentConfiguration!(femmodel, analysis)

   #Call solution sequence to compute new speeds
	println("   extrapolation of ", EnumToString(extvar_num), ":");
   solutionsequence_linear(femmodel,analysis)

	# save
	RequestedOutputsx(femmodel, [extvar_num])
	return nothing
end #}}}
function CreatePVector(analysis::ExtrapolationAnalysis,element::Tria)# {{{
	return nothing
end #}}}
function CreateKMatrix(analysis::ExtrapolationAnalysis,element::Tria)# {{{

	#Internmediaries
	numnodes = 3
	
	#Initialize Element matrix and basis function derivatives
	Ke = ElementMatrix(element.nodes)
	dbasis = Matrix{Float64}(undef,numnodes,2)
	basis  = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	#lsf_slopex_input = GetInput(element, LevelsetfunctionSlopeXEnum)
	#lsf_slopey_input = GetInput(element, LevelsetfunctionSlopeYEnum)

	#Start integrating
	gauss = GaussTria(P1Enum)
	for ig in 1:gauss.numgauss

		Jdet = JacobianDeterminant(xyz_list, gauss)
		NodalFunctionsDerivatives(element, dbasis, xyz_list, gauss)
		NodalFunctions(element, basis, gauss, ig, P1Enum)

		D_scalar = gauss.weights[ig]*Jdet

		for i in 1:numnodes
			for j in 1:numnodes
				Ke.values[i ,j] += D_scalar * (dbasis[j,1]*dbasis[i,1] + dbasis[j,2]*dbasis[i,2])
			end
		end
	end

	return Ke
end #}}}
function GetSolutionFromInputs(analysis::ExtrapolationAnalysis,ug::IssmVector,element::Tria) #{{{
	return nothing
end#}}}
function InputUpdateFromSolution(analysis::ExtrapolationAnalysis,ug::Vector{Float64},element::Tria) #{{{
	extvar_num  = FindParam(IssmEnum, element, ExtrapolationVariableEnum)
	InputUpdateFromSolutionOneDof(element, ug, extvar_num)
end#}}}
function SetConstraintsOnIce(analysis::ExtrapolationAnalysis, element::Tria) #{{{

	# Intermediaries
	numnodes = 3

	# get parameters
	extvar_num  = FindParam(IssmEnum, element, ExtrapolationVariableEnum)

	active_input  = GetInput(element, IceMaskNodeActivationEnum)
	extvar_input  = GetInput(element, extvar_num)

	gauss = GaussTria(P1Enum)
	for ig in 1:gauss.numgauss
		node = element.nodes[ig]
		active = GetInputValue(active_input, gauss, ig)
		if node.active
			if active>0.5
				value = GetInputValue(extvar_input, gauss, ig)
				ApplyConstraint!(node, Int8(1), value)
			else
				DofInFSet!(node, Int8(1))
			end
		end
	end
	return nothing
end#}}}
function UpdateConstraints(analysis::ExtrapolationAnalysis, femmodel::FemModel) #{{{
	for i=1:length(femmodel.elements)
		SetConstraintsOnIce(analysis,femmodel.elements[i])
	end
	return nothing
end#}}}
