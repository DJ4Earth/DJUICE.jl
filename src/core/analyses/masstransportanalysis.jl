#MasstransportAnalysis class definition
struct MasstransportAnalysis <: Analysis#{{{
end #}}}

#Model Processing
function CreateConstraints(analysis::MasstransportAnalysis,constraints::Vector{Constraint},md::model) #{{{

	#load constraints from model
	spcthickness = md.masstransport.spcthickness
	@assert size(spcthickness,1)==md.mesh.numberofvertices
	@assert size(spcthickness,2)==1

	count = 1
	for i in 1:md.mesh.numberofvertices
		if ~isnan(spcthickness[i])
			push!(constraints,Constraint(count,i,1,spcthickness[i]))
			count+=1
		end
	end

	return nothing
end#}}}
function CreateNodes(analysis::MasstransportAnalysis,nodes::Vector{Node},md::model) #{{{

	numdof = 1
	for i in 1:md.mesh.numberofvertices
		push!(nodes,Node(i,i,true,true,numdof,-ones(Int64,numdof), ones(Int64,numdof), -ones(Int64,numdof), zeros(numdof)))
	end

	return nothing
end#}}}
function UpdateElements(analysis::MasstransportAnalysis,elements::Vector{Tria}, inputs::Inputs, md::model) #{{{

	#Provide node indices to element
	for i in 1:md.mesh.numberofelements
		Update(elements[i],inputs,i,md,P1Enum)
	end

	#Add necessary inputs to perform this analysis
	FetchDataToInput(md,inputs,elements,md.geometry.thickness,ThicknessEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.surface,SurfaceEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.base,BaseEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.bed,BedEnum)
	FetchDataToInput(md,inputs,elements,md.basalforcings.groundedice_melting_rate,BasalforcingsGroundediceMeltingRateEnum, 1.0/md.constants.yts)
	FetchDataToInput(md,inputs,elements,md.basalforcings.floatingice_melting_rate,BasalforcingsFloatingiceMeltingRateEnum, 1.0/md.constants.yts)
	FetchDataToInput(md,inputs,elements,md.smb.mass_balance,SmbMassBalanceEnum, 1.0/md.constants.yts)
	FetchDataToInput(md,inputs,elements,md.mask.ice_levelset, MaskIceLevelsetEnum)
	FetchDataToInput(md,inputs,elements,md.mask.ocean_levelset, MaskOceanLevelsetEnum)
	FetchDataToInput(md,inputs,elements,md.initialization.vx,VxEnum, 1.0/md.constants.yts)
	FetchDataToInput(md,inputs,elements,md.initialization.vy,VyEnum, 1.0/md.constants.yts)

	return nothing
end#}}}
function UpdateParameters(analysis::MasstransportAnalysis,parameters::Parameters,md::model) #{{{

	AddParam(parameters, md.masstransport.min_thickness, MasstransportMinThicknessEnum)
	AddParam(parameters, md.masstransport.stabilization, MasstransportStabilizationEnum)

	return nothing
end#}}}

#Finite Element Analysis
function Core(analysis::MasstransportAnalysis,femmodel::FemModel)# {{{

	println("   computing mass transport")
	SetCurrentConfiguration!(femmodel, analysis)

	InputDuplicatex(femmodel, ThicknessEnum, ThicknessOldEnum)
	InputDuplicatex(femmodel, BaseEnum, BaseOldEnum)
	InputDuplicatex(femmodel, SurfaceEnum, SurfaceOldEnum)

	solutionsequence_linear(femmodel,analysis)

	#Save output
	RequestedOutputsx(femmodel, [ThicknessEnum, SurfaceEnum, BaseEnum])

	return nothing
end #}}}
function CreateKMatrix(analysis::MasstransportAnalysis,element::Tria)# {{{
	stabilization = FindParam(Int64, element, MasstransportStabilizationEnum)
	CreateKMatrix(analysis, element, Val(stabilization))::DJUICE.ElementMatrix
end #}}}
function CreateKMatrix(analysis::MasstransportAnalysis,element::Tria, ::Val{stabilization}) where stabilization# {{{

	#Internmediaries
	numnodes = 3

	#Initialize Element matrix and basis function derivatives
	Ke = ElementMatrix(element.nodes)::DJUICE.ElementMatrix
	dbasis = Matrix{Float64}(undef,numnodes,2)
	basis  = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	vx_input      = GetInput(element, VxEnum)
	vy_input      = GetInput(element, VyEnum)
	dt            = FindParam(Float64, element, TimesteppingTimeStepEnum)

	h = CharacteristicLength(element)

	#Start integrating
	gauss = GaussTria(2)
	for ig in 1:gauss.numgauss

		Jdet = JacobianDeterminant(xyz_list, gauss)
		NodalFunctionsDerivatives(element, dbasis, xyz_list, gauss)
		NodalFunctions(element, basis, gauss, ig, P1Enum)

		#Transient term
		for i in 1:numnodes
			for j in 1:numnodes
				Ke.values[i ,j] += gauss.weights[ig]*Jdet*basis[i]*basis[j]
			end
		end

		#Advection term
		vx  = GetInputValue(vx_input, gauss, ig)
		vy  = GetInputValue(vy_input, gauss, ig)
		dvx = GetInputDerivativeValue(vx_input, xyz_list, gauss, ig)
		dvy = GetInputDerivativeValue(vy_input, xyz_list, gauss, ig)
		for i in 1:numnodes
			for j in 1:numnodes
				#\phi_i \phi_j \nabla\cdot v
				Ke.values[i ,j] += dt*gauss.weights[ig]*Jdet*basis[i]*basis[j]*(dvx[1] + dvy[2])
				#\phi_i v\cdot\nabla\phi_j
				Ke.values[i ,j] += dt*gauss.weights[ig]*Jdet*basis[i]*(vx*dbasis[j,1] + vy*dbasis[j,2])
			end
		end

		#Stabilization
		if(stabilization==0)
			#do nothing
		elseif (stabilization==1)
			vx = GetInputAverageValue(vx_input)
			vy = GetInputAverageValue(vy_input)
			D_scalar  = dt*gauss.weights[ig]*Jdet*h*0.5
			# D= | abs(vx),     0 |
			#    | 0,      abs(vy)|
			for i in 1:numnodes 
				for j in 1:numnodes
					Ke.values[i ,j] += D_scalar*(dbasis[i,1]*(abs(vx)*dbasis[j,1]) +
														  dbasis[i,2]*(abs(vy)*dbasis[j,2]))
				end 
			end
		else
			error("Stabilization ",stabilization, " not supported yet")
		end
	end

	return Ke
end #}}}
function CreatePVector(analysis::MasstransportAnalysis,element::Tria)# {{{
	#Internmediaries
	numnodes = 3

	#Initialize Element vectro and basis functions
	pe = ElementVector(element.nodes)
	basis = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	H_input         = GetInput(element, ThicknessEnum)
	gmb_input       = GetInput(element, BasalforcingsGroundediceMeltingRateEnum)
	fmb_input       = GetInput(element, BasalforcingsFloatingiceMeltingRateEnum)
	smb_input       = GetInput(element, SmbMassBalanceEnum)
	olevelset_input = GetInput(element, MaskOceanLevelsetEnum)
	dt            = FindParam(Float64, element, TimesteppingTimeStepEnum)
	stabilization = FindParam(Int64, element, MasstransportStabilizationEnum)

	#How much is actually grounded?
	phi=GetGroundedPortion(element, xyz_list)

	#Start integrating
	gauss = GaussTria(3)
	for ig in 1:gauss.numgauss

		Jdet = JacobianDeterminant(xyz_list, gauss)
		NodalFunctions(element, basis, gauss, ig, P1Enum)

		smb = GetInputValue(smb_input, gauss, ig)
		H   = GetInputValue(H_input, gauss, ig)

		#Only apply melt on fully floating cells
		if(phi<0.00000001)
			mb = GetInputValue(fmb_input, gauss, ig)
		else
			mb = GetInputValue(gmb_input, gauss, ig)
		end

		for i in 1:numnodes
			pe.values[i] += gauss.weights[ig]*Jdet*(H + dt*(smb - mb))*basis[i]
		end
	end

	return pe
end #}}}
function GetSolutionFromInputs(analysis::MasstransportAnalysis,ug::IssmVector,element::Tria) #{{{

	#Get dofs for this finite element
	doflist = GetDofList(element,GsetEnum)
	@assert length(doflist)==3

	#Fetch inputs
	thickness_input = GetInput(element, ThicknessEnum)

	#Loop over each node and enter solution in ug
	count = 0
	gauss=GaussTria(P1Enum)
	for i in 1:gauss.numgauss
		thickness = GetInputValue(thickness_input, gauss, i)

		count += 1
		ug.vector[doflist[count]] = thickness
	end

	#Make sure we reached all the values
	@assert count==length(doflist)

	return nothing
end#}}}
function InputUpdateFromSolution(analysis::MasstransportAnalysis,ug::Vector{Float64},element::Tria) #{{{

	#Get dofs for this finite element
	doflist = GetDofList(element,GsetEnum)

	#Get solution vector for this element
	numdof   = 3
	values = Vector{Float64}(undef,numdof)
	for i in 1:numdof values[i]=ug[doflist[i]] end

	#Get some parameters
	rho_water = FindParam(Float64, element, MaterialsRhoSeawaterEnum)
	rho_ice   = FindParam(Float64, element, MaterialsRhoIceEnum)
	H_min     = FindParam(Float64, element, MasstransportMinThicknessEnum)

	#Now split solution vector into x and y components
	numnodes = 3
	thickness  = Vector{Float64}(undef,numnodes)
	for i in 1:numnodes 
		thickness[i]=values[i] 
		@assert isfinite(thickness[i])

		#Enforce minimum thickness
		if(thickness[i]<H_min)
			thickness[i] = H_min
		end
	end
	AddInput(element, ThicknessEnum,  thickness,  P1Enum)

	#Update bed and surface accordingly
	newthickness = Vector{Float64}(undef,3)
	oldthickness = Vector{Float64}(undef,3)
	oldbase      = Vector{Float64}(undef,3)
	oldsurface   = Vector{Float64}(undef,3)
	phi          = Vector{Float64}(undef,3)
	bed          = Vector{Float64}(undef,3)
	GetInputListOnVertices!(element, newthickness, ThicknessEnum)
	GetInputListOnVertices!(element, oldthickness, ThicknessOldEnum)
	GetInputListOnVertices!(element, oldbase, BaseOldEnum)
	GetInputListOnVertices!(element, oldsurface, SurfaceOldEnum)
	GetInputListOnVertices!(element, phi, MaskOceanLevelsetEnum)
	GetInputListOnVertices!(element, bed, BedEnum)
	sealevel = zeros(3)
	newsurface = Vector{Float64}(undef,3)
	newbase    = Vector{Float64}(undef,3)

	for i in 1:3
		if(phi[i]>0.)
			#this is grounded ice: just add thickness to base.
			newsurface[i] = bed[i]+newthickness[i] #surface = bed + newthickness
			newbase[i]    = bed[i]                 #new base at new bed
		else
			#this is an ice shelf: hydrostatic equilibrium
			newsurface[i] = newthickness[i]*(1-rho_ice/rho_water) + sealevel[i]
			newbase[i]    = newthickness[i]*(-rho_ice/rho_water)  + sealevel[i]
		end
	end

	AddInput(element, SurfaceEnum, newsurface, P1Enum)
	AddInput(element, BaseEnum,    newbase,    P1Enum)

	return nothing
end#}}}
function UpdateConstraints(analysis::MasstransportAnalysis, femmodel::FemModel) #{{{
	SetActiveNodesLSMx(femmodel)
	return nothing
end#}}}
