#LevelsetAnalysis class definition
struct LevelsetAnalysis <: Analysis#{{{
end #}}}

#Model Processing
function CreateConstraints(analysis::LevelsetAnalysis,constraints::Vector{Constraint},md::model) #{{{

	#load constraints from model
	spclevelset = md.levelset.spclevelset

	count = 1
	for i in 1:md.mesh.numberofvertices
		if ~isnan(spclevelset[i])
			push!(constraints,Constraint(count,i,1,spclevelset[i]))
			count+=1
		end
	end

	return nothing
end#}}}
function CreateNodes(analysis::LevelsetAnalysis,nodes::Vector{Node},md::model) #{{{

	numdof = 2
	for i in 1:md.mesh.numberofvertices
		push!(nodes,Node(i,i,true,true,numdof,-ones(Int64,numdof), ones(Int64,numdof), -ones(Int64,numdof), zeros(numdof)))
	end

	return nothing
end#}}}
function UpdateElements(analysis::LevelsetAnalysis,elements::Vector{Tria}, inputs::Inputs, md::model) #{{{

	#Provide node indices to element
	for i in 1:md.mesh.numberofelements
		Update(elements[i],inputs,i,md,P1Enum)
	end

	#Add necessary inputs to perform this analysis
	FetchDataToInput(md,inputs,elements,md.mask.ice_levelset,MaskIceLevelsetEnum)
	FetchDataToInput(md,inputs,elements,md.mask.ocean_levelset,MaskOceanLevelsetEnum)
	FetchDataToInput(md,inputs,elements,md.initialization.vx./md.constants.yts,VxEnum)
	FetchDataToInput(md,inputs,elements,md.initialization.vy./md.constants.yts,VyEnum)

	FetchDataToInput(md,inputs,elements,md.geometry.thickness,ThicknessEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.surface,SurfaceEnum)
	FetchDataToInput(md,inputs,elements,md.geometry.base,BaseEnum)
	FetchDataToInput(md,inputs,elements,md.mask.ice_levelset, MaskIceLevelsetEnum)
	FetchDataToInput(md,inputs,elements,md.mask.ocean_levelset, MaskOceanLevelsetEnum)

	#Get moving front parameters
	if typeof(md.calving) == DefaultCalving
		FetchDataToInput(md,inputs,elements,md.calving.calvingrate,CalvingCalvingrateEnum)
	else
		error("Calving ", typeof(md.calving), " not supported yet")
	end

	#Get frontal melt parameters
	FetchDataToInput(md,inputs,elements,md.frontalforcings.meltingrate, CalvingMeltingrateEnum)
	FetchDataToInput(md,inputs,elements,md.frontalforcings.ablationrate, CalvingAblationrateEnum)

	return nothing
end#}}}
function UpdateParameters(analysis::LevelsetAnalysis,parameters::Parameters,md::model) #{{{
	AddParam(parameters,md.levelset.stabilization,LevelsetStabilizationEnum)
	AddParam(parameters,md.levelset.reinit_frequency,LevelsetReinitFrequencyEnum)
	AddParam(parameters,md.levelset.kill_icebergs,LevelsetKillIcebergsEnum)
	AddParam(parameters,md.levelset.migration_max,MigrationMaxEnum)

	#Deal with Calving
	if typeof(md.calving)==DefaultCalving
		AddParam(parameters, 1, CalvingLawEnum)
	else
		error("Calving ", typeof(md.calving), " not supported yet")
	end

	return nothing
end#}}}

#Finite Element Analysis
function Core(analysis::LevelsetAnalysis,femmodel::FemModel)# {{{

	# moving front
	MovingFrontalVel(femmodel)
	#Activate formulation
	SetCurrentConfiguration!(femmodel, analysis)

	#Call solution sequence to compute new speeds
	println("   call computational core:");
	solutionsequence_linear(femmodel,analysis)

	return nothing
end #}}}
function CreateKMatrix(analysis::LevelsetAnalysis,element::Tria)# {{{

	#Internmediaries
	numnodes = 3
	
	#Initialize Element matrix and basis function derivatives
	Ke = ElementMatrix(element.nodes)
	dbasis = Matrix{Float64}(undef,numnodes,2)
	basis  = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	mf_vx_input      = GetInput(element, MovingFrontalVxEnum)
	mf_vy_input      = GetInput(element, MovingFrontalVyEnum)
	dt            = FindParam(Float64, element, TimesteppingTimeStepEnum)
	stabilization = FindParam(Int64, element, LevelsetStabilizationEnum)
	migration_max  = FindParam(Float64, element, MigrationMaxEnum)

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

      # levelset speed
		wx  = GetInputValue(mf_vx_input, gauss, ig)
      wy  = GetInputValue(mf_vy_input, gauss, ig)
		vel = sqrt(wx^2+wy^2)+1.0e-14

		#rescale by migration_max
		if (vel > migration_max)
			wx = wx/vel*migration_max
			wy = wy/vel*migration_max
		end

      for i in 1:numnodes
         for j in 1:numnodes
            #\phi_i v\cdot\nabla\phi_j
            Ke.values[i ,j] += dt*gauss.weights[ig]*Jdet*basis[i]*(wx*dbasis[j,1] + wy*dbasis[j,2])
         end
      end

		#Stabilization
		if(stabilization==0)
			#do nothing
		elseif (stabilization==1)
			hx, hy, hz = GetElementSizes(element)
			h = sqrt((hx*wx/vel)^2 + (hy*wy/vel)^2)
			kappa = h*vel/2.0
			D  = dt*gauss.weights[ig]*Jdet*kappa
			for i in 1:numnodes
				for j in 1:numnodes
					for k in 1:2
						Ke.values[i ,j] += D*dbasis[i,k]*dbasis[j,k]
					end
				end 
			end
		else
			error("Stabilization ",stabilization, " not supported yet")
		end
	end

	return Ke
end #}}}
function CreatePVector(analysis::LevelsetAnalysis,element::Tria)# {{{

	#Internmediaries
	numnodes = 3

	#Initialize Element vectro and basis functions
	pe = ElementVector(element.nodes)
	basis = Vector{Float64}(undef,numnodes)

	#Retrieve all inputs and parameters
	xyz_list = GetVerticesCoordinates(element.vertices)
	levelset_input      = GetInput(element, MaskIceLevelsetEnum)
	mf_vx_input      = GetInput(element, MovingFrontalVxEnum)
	mf_vy_input      = GetInput(element, MovingFrontalVyEnum)
	dt            = FindParam(Float64, element, TimesteppingTimeStepEnum)
	stabilization = FindParam(Int64, element, LevelsetStabilizationEnum)
	migration_max  = FindParam(Float64, element, MigrationMaxEnum)

	h = CharacteristicLength(element)

	#Start integrating
	gauss = GaussTria(2)
	for ig in 1:gauss.numgauss

		Jdet = JacobianDeterminant(xyz_list, gauss)
		#Get nodal basis
		NodalFunctions(element, basis, gauss, ig, P1Enum)
		#Old function value
		lsf = GetInputValue(levelset_input, gauss, ig)
		for i in 1:numnodes
         pe.values[i] += gauss.weights[ig]*Jdet*lsf*basis[i]
		end
		#TODO: add stab=5
	end

	return pe
end #}}}
function GetSolutionFromInputs(analysis::LevelsetAnalysis,ug::IssmVector,element::Tria) #{{{

	#Get dofs for this finite element
	doflist = GetDofList(element,GsetEnum)
	@assert length(doflist)==3

	#Fetch inputs
	mask_input = GetInput(element, MaskIceLevelsetEnum)

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
function InputUpdateFromSolution(analysis::LevelsetAnalysis,ug::Vector{Float64},element::Tria) #{{{
	InputUpdateFromSolutionOneDof(element, ug, MaskIceLevelsetEnum)
end#}}}
function UpdateConstraints(analysis::LevelsetAnalysis, femmodel::FemModel) #{{{
	#Default do nothing
	return nothing
end#}}}

# Moving front
function MovingFrontalVel(femmodel::FemModel)# {{{
	for i in 1:length(femmodel.elements)
		MovingFrontalVelocity(femmodel.elements[i])
	end
	return nothing
end #}}}
