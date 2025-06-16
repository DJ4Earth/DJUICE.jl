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
   extvar_num  = FindParam(Float64, femmodel.parameters,ExtrapolationVariableEnum)

	SetCurrentConfiguration!(femmodel, analysis)

   #Call solution sequence to compute new speeds
	println("   extrapolation of", EnumToString(extvar), ":");
   solutionsequence_linear(femmodel,analysis)

	# save
   RequestedOutputsx(femmodel, extvar_num)
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
	lsf_slopex_input = GetInput(element, LevelsetfunctionSlopeXEnum)
	lsf_slopey_input = GetInput(element, LevelsetfunctionSlopeYEnum)

	#Start integrating
	gauss = GaussTria(2)
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
	extvar_num  = FindParam(Float64, femmodel.parameters,ExtrapolationVariableEnum)
	InputUpdateFromSolutionOneDof(element, ug, extvar_num)
end#}}}
function UpdateConstraints(analysis::ExtrapolationAnalysis, femmodel::FemModel) #{{{
	#Default do nothing
	return nothing
end#}}}
