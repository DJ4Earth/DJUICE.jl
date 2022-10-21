#StressbalanceAnalysis class definition
struct StressbalanceAnalysis <: Analysis#{{{
end #}}}

#Model Processing
function CreateConstraints(analysis::StressbalanceAnalysis,constraints::Vector{Constraint},md::model) #{{{

	#load constraints from model
	spcvx = md.stressbalance.spcvx
	spcvy = md.stressbalance.spcvy

	count = 1
	for i in 1:md.mesh.numberofvertices
		if ~isnan(spcvx[i])
			push!(constraints,Constraint(count,i,1,spcvx[i]/md.constants.yts))
			count+=1
		end
		if ~isnan(spcvy[i])
			push!(constraints,Constraint(count,i,2,spcvy[i]/md.constants.yts))
			count+=1
		end
	end

	return nothing
end#}}}
function CreateNodes(analysis::StressbalanceAnalysis,nodes::Vector{Node},md::model) #{{{

	numdof = 2
	for i in 1:md.mesh.numberofvertices
		push!(nodes,Node(i,i,true,true,numdof,-ones(Int64,numdof), ones(Int64,numdof), -ones(Int64,numdof), zeros(numdof)))
	end

	return nothing
end#}}}
function UpdateElements(analysis::StressbalanceAnalysis,elements::Vector{Tria}, inputs::Inputs, md::model) #{{{

	#Provide node indices to element
	for i in 1:md.mesh.numberofelements
		Update(elements[i],inputs,i,md,P1Enum)
	end

	#Add necessary inputs to perform this analysis
	FetchDataToInput(md,inputs,elements,md.materials.rheology_B,MaterialsRheologyBEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.thickness,ThicknessEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.surface,SurfaceEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.base,BaseEnum)
	FetchDataToInput(md,inputs,elements,md.initialization.vx./md.constants.yts,VxEnum)
	FetchDataToInput(md,inputs,elements,md.initialization.vy./md.constants.yts,VyEnum)
	FetchDataToInput(md,inputs,elements,md.mask.ice_levelset, MaskIceLevelsetEnum)
	FetchDataToInput(md,inputs,elements,md.mask.ocean_levelset, MaskOceanLevelsetEnum)

	#Deal with friction
	if typeof(md.friction) == BuddFriction
		FetchDataToInput(md,inputs,elements,md.friction.coefficient,FrictionCoefficientEnum)
		FetchDataToInput(md,inputs,elements,md.friction.p,FrictionPEnum)
		FetchDataToInput(md,inputs,elements,md.friction.q,FrictionQEnum)
	elseif typeof(md.friction) == WeertmanFriction
		FetchDataToInput(md,inputs,elements,md.friction.C,FrictionCEnum)
		FetchDataToInput(md,inputs,elements,md.friction.m,FrictionMEnum)
	elseif typeof(md.friction) == DNNFriction
		FetchDataToInput(md,inputs,elements,md.friction.coefficient,FrictionCoefficientEnum)
		FetchDataToInput(md,inputs,elements,md.geometry.ssx,SurfaceSlopeXEnum)
		FetchDataToInput(md,inputs,elements,md.geometry.ssy,SurfaceSlopeYEnum)
		FetchDataToInput(md,inputs,elements,md.geometry.bsx,BedSlopeXEnum)
		FetchDataToInput(md,inputs,elements,md.geometry.bsy,BedSlopeYEnum)
	else
		error("Friction ", typeof(md.friction), " not supported yet")
	end

	return nothing
end#}}}
function UpdateParameters(analysis::StressbalanceAnalysis,parameters::Parameters,md::model) #{{{
	AddParam(parameters,md.stressbalance.restol,StressbalanceRestolEnum)
	AddParam(parameters,md.stressbalance.reltol,StressbalanceReltolEnum)
	AddParam(parameters,md.stressbalance.abstol,StressbalanceAbstolEnum)
	AddParam(parameters,md.stressbalance.maxiter,StressbalanceMaxiterEnum)

	#Deal with friction
	if typeof(md.friction)==BuddFriction
		AddParam(parameters, 1, FrictionLawEnum)
	elseif typeof(md.friction)==WeertmanFriction
		AddParam(parameters, 2, FrictionLawEnum)
	elseif typeof(md.friction)==DNNFriction
		AddParam(parameters, 10, FrictionLawEnum)
		AddParam(parameters, md.friction.dnnChain, FrictionDNNChainEnum)
		AddParam(parameters, md.friction.dtx, FrictionDNNdtxEnum)
		AddParam(parameters, md.friction.dty, FrictionDNNdtyEnum)
	else
		error("Friction ", typeof(md.friction), " not supported yet")
	end

	return nothing
end#}}}

#Finite Element Analysis
function Core(analysis::StressbalanceAnalysis,femmodel::FemModel)# {{{

	#Set current analysis to Stressnalance
	SetCurrentConfiguration!(femmodel, analysis)

	#Fetch parameters relevant to solution sequence
	maxiter = FindParam(Int64, femmodel.parameters,StressbalanceMaxiterEnum)
	restol  = FindParam(Float64, femmodel.parameters,StressbalanceRestolEnum)
	reltol  = FindParam(Float64, femmodel.parameters,StressbalanceReltolEnum)
	abstol  = FindParam(Float64, femmodel.parameters,StressbalanceAbstolEnum)

	#Call solution sequence to compute new speeds
	println("   computing stress balance");
	solutionsequence_nonlinear(femmodel,analysis,maxiter,restol,reltol,abstol)

	#Save output
	RequestedOutputsx(femmodel, [VxEnum,VyEnum,VelEnum])

	return nothing
end #}}}
function CreateKMatrix(analysis::StressbalanceAnalysis,element::Tria)# {{{

	if(!IsIceInElement(element)) return end

	#Internmediaries
	numnodes = 3
	
	#Initialize Element matrix and basis function derivatives
	Ke = ElementMatrix(element.nodes)
	dbasis = Matrix{Float64}(undef,numnodes,2)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	H_input  = GetInput(element, ThicknessEnum)

	#Prepare material object
	material = Matice(element)
	
	#Start integrating
	gauss = GaussTria(2)
	for ig in 1:gauss.numgauss

		Jdet = JacobianDeterminant(xyz_list, gauss)
		NodalFunctionsDerivatives(element,dbasis,xyz_list,gauss)

		H  = GetInputValue(H_input, gauss, ig)
		mu = ViscositySSA(material, xyz_list, gauss, ig)

		for i in 1:numnodes
			for j in 1:numnodes
				Ke.values[2*i-1,2*j-1] += gauss.weights[ig]*Jdet*mu*H*(4*dbasis[j,1]*dbasis[i,1] + dbasis[j,2]*dbasis[i,2])
				Ke.values[2*i-1,2*j  ] += gauss.weights[ig]*Jdet*mu*H*(2*dbasis[j,2]*dbasis[i,1] + dbasis[j,1]*dbasis[i,2])
				Ke.values[2*i  ,2*j-1] += gauss.weights[ig]*Jdet*mu*H*(2*dbasis[j,1]*dbasis[i,2] + dbasis[j,2]*dbasis[i,1])
				Ke.values[2*i  ,2*j  ] += gauss.weights[ig]*Jdet*mu*H*(4*dbasis[j,2]*dbasis[i,2] + dbasis[j,1]*dbasis[i,1])
			end
		end
	end

	#Add basal friction
	phi=GetGroundedPortion(element, xyz_list)

	if(phi>0)
		basis = Vector{Float64}(undef,numnodes)
		friction = CoreFriction(element)

		#Start integrating
		gauss = GaussTria(2)
		for ig in 1:gauss.numgauss

			Jdet = JacobianDeterminant(xyz_list, gauss)
			NodalFunctions(element, basis, gauss, ig, P1Enum)

			alpha2 = Alpha2(friction, gauss, ig)

			for i in 1:numnodes
				for j in 1:numnodes
					Ke.values[2*i-1,2*j-1] += gauss.weights[ig]*Jdet*phi*alpha2*basis[i]*basis[j]
					Ke.values[2*i  ,2*j  ] += gauss.weights[ig]*Jdet*phi*alpha2*basis[i]*basis[j]
				end
			end
		end
	end

	return Ke
end #}}}
function CreatePVector(analysis::StressbalanceAnalysis,element::Tria)# {{{

	if(!IsIceInElement(element)) return end

	#Internmediaries
	numnodes = 3

	#Initialize Element vectro and basis functions
	pe = ElementVector(element.nodes)
	basis = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	H_input  = GetInput(element, ThicknessEnum)
	s_input  = GetInput(element, SurfaceEnum)
	rho_ice  = FindParam(Float64, element, MaterialsRhoIceEnum)
	g        = FindParam(Float64, element, ConstantsGEnum)

	#Start integrating
	gauss = GaussTria(2)
	for ig in 1:gauss.numgauss

		Jdet = JacobianDeterminant(xyz_list, gauss)
		NodalFunctions(element, basis, gauss, ig, P1Enum)

		H  = GetInputValue(H_input, gauss, ig)
		ds = GetInputDerivativeValue(s_input, xyz_list, gauss, ig)

		for i in 1:numnodes
			pe.values[2*i-1] += -gauss.weights[ig]*Jdet*rho_ice*g*H*ds[1]*basis[i]
			pe.values[2*i  ] += -gauss.weights[ig]*Jdet*rho_ice*g*H*ds[2]*basis[i]
		end
	end

	if(IsIcefront(element))

		#Get additional parameters and inputs
		b_input   = GetInput(element, BaseEnum)
		rho_water = FindParam(Float64, element, MaterialsRhoSeawaterEnum)

		#Get normal and ice front coordinates
		xyz_list_front = Matrix{Float64}(undef,2,3)
		GetIcefrontCoordinates!(element, xyz_list_front, xyz_list, MaskIceLevelsetEnum)
		nx, ny = NormalSection(element, xyz_list_front)

		gauss = GaussTria(element, xyz_list, xyz_list_front, 3)
		for ig in 1:gauss.numgauss

			Jdet = JacobianDeterminantSurface(xyz_list_front, gauss)
			NodalFunctions(element, basis, gauss, ig, P1Enum)

			H  = GetInputValue(H_input, gauss, ig)
			b  = GetInputValue(b_input, gauss, ig)
			sl = 0

			term = 0.5*g*rho_ice*H^2 + 0.5*g*rho_water*(min(0, H+b-sl)^2 - min(0, b-sl)^2)

			for i in 1:numnodes
				pe.values[2*i-1] += gauss.weights[ig]*Jdet*term*nx*basis[i]
				pe.values[2*i  ] += gauss.weights[ig]*Jdet*term*ny*basis[i]
			end
		end
	end

	return pe
end #}}}
function GetSolutionFromInputs(analysis::StressbalanceAnalysis,ug::IssmVector,element::Tria) #{{{

	#Get dofs for this finite element
	doflist = GetDofList(element,GsetEnum)
	@assert length(doflist)==6

	#Fetch inputs
	vx_input = GetInput(element, VxEnum)
	vy_input = GetInput(element, VyEnum)

	#Loop over each node and enter solution in ug
	count = 0
	gauss=GaussTria(P1Enum)
	for i in 1:gauss.numgauss
		vx = GetInputValue(vx_input, gauss, i)
		vy = GetInputValue(vy_input, gauss, i)

		count += 1
		ug.vector[doflist[count]] = vx
		count += 1
		ug.vector[doflist[count]] = vy
	end

	#Make sure we reached all the values
	@assert count==length(doflist)

	return nothing
end#}}}
function InputUpdateFromSolution(analysis::StressbalanceAnalysis,ug::Vector{Float64},element::Tria) #{{{

	#Get dofs for this finite element
	doflist = GetDofList(element,GsetEnum)

	#Get solution vector for this element
	numdof   = 3*2
	values = Vector{Float64}(undef,numdof)
	for i in 1:numdof values[i]=ug[doflist[i]] end

	#Now split solution vector into x and y components
	numnodes = 3
	vx  = Vector{Float64}(undef,numnodes)
	vy  = Vector{Float64}(undef,numnodes)
	vel = Vector{Float64}(undef,numnodes)
	for i in 1:numnodes 
		vx[i]=values[2*i-1] 
		vy[i]=values[2*i] 
		@assert isfinite(vx[i])
		@assert isfinite(vy[i])

		vel[i] =sqrt(vx[i]^2 + vy[i]^2)
	end

	AddInput(element, VxEnum,  vx,  P1Enum)
	AddInput(element, VyEnum,  vy,  P1Enum)
	AddInput(element, VelEnum, vel, P1Enum)

	return nothing
end#}}}
function UpdateConstraints(analysis::StressbalanceAnalysis, femmodel::FemModel) #{{{
	SetActiveNodesLSMx(femmodel)
	return nothing
end#}}}
