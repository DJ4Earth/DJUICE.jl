#Tria class definition
mutable struct Tria #{{{
	sid::Int64
	pid::Int64

	#vertexids::Int64[3]
	#vertices::Vertex[3]
	vertexids::Vector{Int64}
	vertices::Vector{Vertex}

	nodes::Vector{Node}
	nodes_list::Vector{Vector{Node}}
	nodes_ids_list::Vector{Vector{Int64}}

	parameters::Parameters
	inputs::Inputs
end# }}}
function Tria(sid::Int64, pid::Int64, vertexids::Vector{Int64}) #{{{

	#This is the default constructor, at this point we don't have much information
	tempparams   = Parameters()
	tempinputs   = Inputs(-1,-1,Dict{IssmEnum,Input}())
	return Tria(sid, pid,
					vertexids, Vector{Vertex}(undef,3),
					Vector{Node}(undef,0), Vector{Vector{Node}}(undef,0), Vector{Vector{Int64}}(undef,0),
					tempparams, tempinputs)

end #}}}

#Element functions
function InputCreate(element::Tria,inputs::Inputs,data::Vector{Float64},enum::IssmEnum) #{{{
	if size(data,1)==inputs.numberofelements
		SetTriaInput(inputs,enum,P0Enum,element.sid,data[element.sid])
	elseif size(data,1)==inputs.numberofvertices
		SetTriaInput(inputs,enum,P1Enum,element.vertexids,data[element.vertexids])
	else
		error("size ",size(data,1)," not supported for ", enum);
	end

	return nothing
end #}}}
function AddInput(element::Tria,inputenum::IssmEnum,data::Vector{Float64},interpolation::IssmEnum) #{{{
	if interpolation==P1Enum
		@assert length(data)==3
		SetTriaInput(element.inputs,inputenum,P1Enum,element.vertexids,data)
	else
		error("interpolation ", interpolation, " not supported yet");
	end

	return nothing
end #}}}
function InputUpdateFromVector(element::Tria, vector::Vector{Float64}, enum::IssmEnum, layout::IssmEnum) #{{{

	lidlist = element.vertexids
	data = Vector{Float64}(undef, 3)

	if(layout==VertexSIdEnum)
		for i in 1:3
			data[i] = vector[element.vertices[i].sid]
			@assert isfinite(data[i])
		end
		SetTriaInput(element.inputs, enum, P1Enum, lidlist, data)
	else
		error("layout ", layout, " not supported yet");
	end

	return nothing
end #}}}
function Update(element::Tria, inputs::Inputs, index::Int64, md::model, finiteelement::IssmEnum) #{{{

	if finiteelement==P1Enum
		numnodes = 3
		nodeids    = Vector{Int64}(undef,numnodes)
		nodeids[1] = md.mesh.elements[index,1]
		nodeids[2] = md.mesh.elements[index,2]
		nodeids[3] = md.mesh.elements[index,3]

		push!(element.nodes_ids_list, nodeids)
		push!(element.nodes_list, Vector{Node}(undef, numnodes))
	else
		error("not supported yet")
	end

	return nothing
end #}}}
function Configure(element::Tria,nodes::Vector{Node},vertices::Vector{Vertex},parameters::Parameters,inputs::Inputs,index::Int64) # {{{

   #Configure vertices
   for i in 1:3
		element.vertices[i] = vertices[element.vertexids[i]]
   end

	#Configure nodes (assuming P1 finite elements)
	nodes_list     = element.nodes_list[index]
	nodes_ids_list = element.nodes_ids_list[index]
	for i in 1:3
		nodes_list[i] = nodes[nodes_ids_list[i]]
	end

	#Point to real datasets
	element.nodes      = element.nodes_list[index]
	element.parameters = parameters
	element.inputs     = inputs

	return nothing
end # }}}
function GetDofList(element::Tria,setenum::IssmEnum) # {{{

	#Define number of nodes
	numnodes = 3

	#Determine size of doflist
	numdofs = 0
	for i in 1:numnodes
		numdofs += GetNumberOfDofs(element.nodes[i],GsetEnum)
	end

	#Allocate doflist vector
	doflist = Vector{Int64}(undef,numdofs)

	#enter dofs in doflist vector
	count = 0
	for i in 1:numnodes
		count = GetDofList(element.nodes[i],doflist,count,GsetEnum)
	end

	return doflist
end # }}}
function GetInput(element::Tria,enum::IssmEnum) # {{{

	input = GetInput(element.inputs,enum)
	InputServe!(element,input)
	return input

end # }}}
function GetInputListOnNodes!(element::Tria, vector::Vector{Float64}, enum::IssmEnum) # {{{

	#Get Input first 
	input = GetInput(element, enum)

	#Get value at each vertex (i.e. P1 Nodes)
	gauss=GaussTria(P1Enum)
	for i in 1:gauss.numgauss
		vector[i] = GetInputValue(input, gauss, i)
	end

	return nothing
end # }}}
function GetInputValue(element::Tria, vector::Vector{Float64}, gauss::GaussTria, ig::Int64, enum::IssmEnum) # {{{

	# Allocate basis functions
	basis = Vector{Float64}(undef, 3)
	# Fetch number of nodes
	numnodes = NumberofNodesTria(enum)
	@assert numnodes <= 3
	# Get basis functions at this point
	NodalFunctions(element, basis, gauss, ig, enum)
	
	# Calculate parameter for this Gauss point
	value::Float64 = 0.0
	for i = 1:3
		value += basis[i] * vector[i]
	end

	return value
end # }}}
function GetInputListOnVertices!(element::Tria, vector::Vector{Float64}, enum::IssmEnum) # {{{

	#Get Input first 
	input = GetInput(element, enum)

	#Get value at each vertex (i.e. P1 Nodes)
	gauss=GaussTria(P1Enum)
	for i in 1:gauss.numgauss
		vector[i] = GetInputValue(input, gauss, i)
	end

	return nothing
end # }}}
function FindParam(::Type{T}, element::Tria, enum::IssmEnum) where T # {{{

	return FindParam(T, element.parameters, enum)

end # }}}
function InputServe!(element::Tria,input::Input) # {{{

	if input.interp==P0Enum
		input.element_values[1] = input.values[element.sid]
	elseif input.interp==P1Enum
		for i in 1:3
			input.element_values[i] = input.values[element.vertices[i].sid]
		end
	else
		error("interpolation ",input.interp," not supported yet")
	end

	return nothing
end # }}}
function GetGroundedPortion(element::Tria, xyz_list::Matrix{Float64}) #{{{

	level = Vector{Float64}(undef,3)
	GetInputListOnVertices!(element, level, MaskOceanLevelsetEnum)

	#Be sure that values are not zero
	epsilon = 1.e-15
	for i in 1:3
		if(level[i]==0.) level[i]=level[i]+epsilon end
	end

	if level[1]>0 && level[2]>0 && level[3]>0
      #Completely grounded
		phi = 1.0
   elseif level[1]<0 && level[2]<0 && level[3]<0
      #Completely floating
      phi = 0.0
   else
		#Partially floating,
		if(level[1]*level[2]>0) #Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2
			s1=level[3]/(level[3]-level[2]);
			s2=level[3]/(level[3]-level[1]);
		elseif(level[2]*level[3]>0) #Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2
			s1=level[1]/(level[1]-level[2]);
			s2=level[1]/(level[1]-level[3]);
		elseif(level[1]*level[3]>0) #Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2
			s1=level[2]/(level[2]-level[1]);
			s2=level[2]/(level[2]-level[3]);
		else
			error("not supposed to be here...")
		end

		if(level[1]*level[2]*level[3]>0)
			phi = s1*s2
		else
			phi = (1-s1*s2)
		end
	end

	return phi
end#}}}
function InputUpdateFromSolutionOneDof(element::Tria, ug::Vector{Float64},enum::IssmEnum) #{{{
	#Get dofs for this finite element
   doflist = GetDofList(element,GsetEnum)

   #Get solution vector for this element
   numdof   = 3
   values = Vector{Float64}(undef,numdof)
   for i in 1:numdof values[i]=ug[doflist[i]] end

   #Add back to Mask
   AddInput(element, enum, value, P1Enum)

   return nothing

end#}}}
function IsIcefront(element::Tria) #{{{

	level = Vector{Float64}(undef,3)
	GetInputListOnVertices!(element, level, MaskIceLevelsetEnum)

	nbice = 0
	for i in 1:3
		if(level[i]<0.) nbice+=1 end
	end

	if(nbice==1)
		return true
	else
		return false
	end
end#}}}
function IsIceInElement(element::Tria) #{{{
	#We consider that an element has ice if at least one of its nodes has a negative level set

	input=GetInput(element, MaskIceLevelsetEnum)

	if GetInputMin(input)<0
		return true
	else
		return false
	end

end#}}}
function GetIcefrontCoordinates!(element::Tria, xyz_front::Matrix{Float64}, xyz_list::Matrix{Float64}, levelsetenum::IssmEnum) #{{{

	#Intermediaries
	level        = Vector{Float64}(undef,3)
	indicesfront = Vector{Int64}(undef,3)

	#Recover value of levelset for all vertices
	GetInputListOnVertices!(element, level, levelsetenum)

	#Get nodes where there is no ice
	num_frontnodes = 0
	for i in 1:3
		if(level[i]>=0.)
			num_frontnodes += 1
			indicesfront[num_frontnodes] = i
		end
	end
	@assert num_frontnodes==2

	#Arrange order of frontnodes such that they are oriented counterclockwise
	NUMVERTICES = 3
	if((NUMVERTICES+indicesfront[1]-indicesfront[2])%NUMVERTICES != NUMVERTICES-1)
		index=indicesfront[1]
		indicesfront[1]=indicesfront[2]
		indicesfront[2]=index
	end

	#Return nodes
	xyz_front[1,:]=xyz_list[indicesfront[1],:]
	xyz_front[2,:]=xyz_list[indicesfront[2],:]
	return nothing
end#}}}
function GetArea(element::Tria)#{{{

	#Get xyz list
	xyz_list = GetVerticesCoordinates(element.vertices)
	x1 = xyz_list[1,1]; y1 = xyz_list[1,2]
	x2 = xyz_list[2,1]; y2 = xyz_list[2,2]
	x3 = xyz_list[3,1]; y3 = xyz_list[3,2]

	@assert x2*y3 - y2*x3 + x1*y2 - y1*x2 + x3*y1 - y3*x1>0
	return (x2*y3 - y2*x3 + x1*y2 - y1*x2 + x3*y1 - y3*x1)/2
end#}}}
function GetElementSizes(element::Tria)#{{{

	#Get xyz list
	xyz_list = GetVerticesCoordinates(element.vertices)
	x1 = xyz_list[1,1]; y1 = xyz_list[1,2]
	x2 = xyz_list[2,1]; y2 = xyz_list[2,2]
	x3 = xyz_list[3,1]; y3 = xyz_list[3,2]

   xmin=xyz_list[1,1]; xmax=xyz_list[1,1];
   ymin=xyz_list[1,2]; ymax=xyz_list[1,2];

	for i in [2,3]
		if(xyz_list[i,1]<xmin) xmin=xyz_list[i,1] end
      if(xyz_list[i,1]>xmax) xmax=xyz_list[i,1] end
      if(xyz_list[i,2]<ymin) ymin=xyz_list[i,2] end
      if(xyz_list[i,2]>ymax) ymax=xyz_list[i,2] end
	end

   hx = xmax-xmin
   hy = ymax-ymin
   hz = 0.

	return hx, hy, hz
end#}}}
function NormalSection(element::Tria, xyz_front::Matrix{Float64}) #{{{

	#Build output pointing vector
	nx =  xyz_front[2,2] - xyz_front[1,2]
	ny = -xyz_front[2,1] + xyz_front[1,1]

	#normalize
	norm = sqrt(nx^2 + ny^2)
	nx = nx/norm
	ny = ny/norm

	return nx, ny
end#}}}
function CharacteristicLength(element::Tria) #{{{

	return sqrt(2*GetArea(element))
end#}}}
function MigrateGroundingLine(element::Tria) #{{{

	h = Vector{Float64}(undef,3)
	s = Vector{Float64}(undef,3)
	b = Vector{Float64}(undef,3)
	r = Vector{Float64}(undef,3)
	phi = Vector{Float64}(undef,3)
	sl = zeros(3)
	GetInputListOnVertices!(element, h, ThicknessEnum)
	GetInputListOnVertices!(element, s, SurfaceEnum)
	GetInputListOnVertices!(element, b, BaseEnum)
	GetInputListOnVertices!(element, r, BedEnum)
	#GetInputListOnVertices(element, sl, SealevelEnum)
	GetInputListOnVertices!(element, phi, MaskOceanLevelsetEnum)


	rho_water   = FindParam(Float64, element, MaterialsRhoSeawaterEnum)
	rho_ice     = FindParam(Float64, element, MaterialsRhoIceEnum)
	density     = rho_ice/rho_water

	for i in 1:3

		if(phi[i]<=0)
			#reground if base is below bed
			if(b[i]<=r[i])
				b[i] = r[i]
				s[i] = b[i]+h[i]
			end
		else
			bed_hydro=-density*h[i]+sl[i];
			if (bed_hydro>r[i])
				#Unground only if the element is connected to the ice shelf
				s[i] = (1-density)*h[i]+sl[i]
				b[i] = -density*h[i]+sl[i]
			end
		end

		#recalculate phi
		phi[i]=h[i]+(r[i]-sl[i])/density
	end

	#Update inputs
	AddInput(element,MaskOceanLevelsetEnum,phi,P1Enum)
	AddInput(element,SurfaceEnum,s,P1Enum)
	AddInput(element,BaseEnum,b,P1Enum)

	return nothing
end#}}}
function GetXcoord(element::Tria, xyz_list::Matrix{Float64}, gauss::GaussTria, ig::Int64) #{{{

	# create a list of x
	x_list = Vector{Float64}(undef,3)

	for i in 1:3
		x_list[i] = xyz_list[i,1]
	end

   # Get value at gauss point
   return GetInputValue(element, x_list, gauss, ig, P1Enum);
end#}}}
function GetYcoord(element::Tria, xyz_list::Matrix{Float64}, gauss::GaussTria, ig::Int64) #{{{

	# create a list of y
	y_list = Vector{Float64}(undef,3)

	for i in 1:3
		y_list[i] = xyz_list[i,2]
	end

   # Get value at gauss point
   return GetInputValue(element, y_list, gauss, ig, P1Enum);
end#}}}

#Finite Element stuff
function JacobianDeterminant(xyz_list::Matrix{Float64}, gauss::GaussTria) #{{{

	#Get Jacobian Matrix
	J = Jacobian(xyz_list)

	#Get its determinant
	Jdet = Matrix2x2Determinant(J)

	#check and return
	if(Jdet<0) error("negative Jacobian Determinant") end
	return Jdet

end#}}}
function JacobianDeterminantSurface(xyz_list::Matrix{Float64}, gauss::GaussTria) #{{{

	x1 = xyz_list[1,1]; y1 = xyz_list[1,2]
	x2 = xyz_list[2,1]; y2 = xyz_list[2,2]
	Jdet = .5*sqrt((x2-x1)^2 + (y2-y1)^2)

	#check and return
	if(Jdet<0) error("negative Jacobian Determinant") end
	return Jdet

end#}}}
function Jacobian(xyz_list::Matrix{Float64}) #{{{

	J = Matrix{Float64}(undef,2,2)

	x1 = xyz_list[1,1]
	y1 = xyz_list[1,2]
	x2 = xyz_list[2,1]
	y2 = xyz_list[2,2]
	x3 = xyz_list[3,1]
	y3 = xyz_list[3,2]

	J[1,1] = .5*(x2-x1)
	J[1,2] = .5*(y2-y1)
	J[2,1] = sqrt(3)/6*(2*x3 -x1 -x2)
	J[2,2] = sqrt(3)/6*(2*y3 -y1 -y2)

	return J
end#}}}
function JacobianInvert(xyz_list::Matrix{Float64}, gauss::GaussTria) #{{{

	#Get Jacobian matrix
	J = Jacobian(xyz_list)

	#Get its determinant
	Jinv = Matrix2x2Invert(J)

	return Jinv
end#}}}
function NodalFunctions(element::Tria,basis::Vector{Float64}, gauss::GaussTria, ig::Int64, finiteelement::IssmEnum) #{{{

	if(finiteelement==P0Enum)
		#Nodal function 1
		basis[1]= 1
	elseif(finiteelement==P1Enum)
		basis[1] = gauss.coords1[ig]
		basis[2] = gauss.coords2[ig]
		basis[3] = gauss.coords3[ig]
	else
		error("Element type ",finiteelement," not supported yet")
	end

	return nothing
end#}}}
function NodalFunctionsDerivatives(element::Tria,dbasis::Matrix{Float64},xyz_list::Matrix{Float64}, gauss::GaussTria) #{{{

	#Get nodal function derivatives in reference element
	dbasis_ref = Matrix{Float64}(undef,3,2)
	NodalFunctionsDerivativesReferenceTria(dbasis_ref,gauss,P1Enum)

	#Get invert of the Jacobian
	Jinv = JacobianInvert(xyz_list,gauss)

	#Build dbasis:
	#[ dNi/dx ] = Jinv * [dNhat_i/dr]
	#[ dNi/dy ] =        [dNhat_i/ds]
	for i in 1:3
		dbasis[i,1] = Jinv[1,1]*dbasis_ref[i,1]+Jinv[1,2]*dbasis_ref[i,2]
		dbasis[i,2] = Jinv[2,1]*dbasis_ref[i,1]+Jinv[2,2]*dbasis_ref[i,2]
	end

	return nothing
end#}}}
function NodalFunctionsDerivativesReferenceTria(dbasis::Matrix{Float64}, gauss::GaussTria, finiteelement::IssmEnum) #{{{

	if(finiteelement==P0Enum)
		#Nodal function 1
		dbasis[1,1]= 0
		dbasis[1,2]= 0

	elseif(finiteelement==P1Enum)
		#Nodal function 1
		dbasis[1,1]= -.5
		dbasis[1,2]= -sqrt(3)/6
		#Nodal function 2
		dbasis[2,1]= .5
		dbasis[2,2]= -sqrt(3)/6
		#Nodal function 3
		dbasis[3,1]= 0
		dbasis[3,2]= sqrt(3)/3
	else
		error("Element type ",finiteelement," not supported yet")
	end

	return nothing
end#}}}
function NumberofNodesTria(finiteelement) #{{{

	if    (finiteelement==P0Enum) return 0
	elseif(finiteelement==P1Enum) return 3
	else
		error("Element type ",finiteelement," not supported yet")
	end

	return nothing
end#}}}
function GaussTria(element::Tria, xyz_list::Matrix{Float64}, xyz_list_front::Matrix{Float64}, order::Int64) #{{{

	area_coordinates = Matrix{Float64}(undef,2,3)
	GetAreaCoordinates!(element, area_coordinates, xyz_list_front, xyz_list)

	return GaussTria(area_coordinates, order)
end# }}}
function GetAreaCoordinates!(element::Tria, area_coordinates::Matrix{Float64}, xyz_zero::Matrix{Float64}, xyz_list::Matrix{Float64})#{{{

	numpoints = size(area_coordinates,1)
	area = GetArea(element)

	#Copy original xyz_list
	xyz_bis=copy(xyz_list)
	for i in 1:numpoints
		for j in 1:3

			#Change appropriate line
			xyz_bis[j,:] = xyz_zero[i,:]

			#Compute area fraction
			area_portion=abs(xyz_bis[2,1]*xyz_bis[3,2] - xyz_bis[2,2]*xyz_bis[3,1] + xyz_bis[1,1]*xyz_bis[2,2] - xyz_bis[1,2]*xyz_bis[2,1] + xyz_bis[3,1]*xyz_bis[1,2] - xyz_bis[3,2]*xyz_bis[1,1])/2
			area_coordinates[i,j] = area_portion/area

			#reinitialize xyz_list
			xyz_bis[j,:] = xyz_list[j,:]
		end
	end

	return nothing
end #}}}
