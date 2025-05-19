function setmask(md::model,floatingicename::String,groundedicename::String) #{{{
	"""
	SETMASK - establish boundaries between grounded and floating ice.

	By default, ice is considered grounded. The contour floatingicename defines nodes 
	for which ice is floating. The contour groundedicename defines nodes inside an floatingice, 
	that are grounded (ie: ice rises, islands, etc ...)
	All input files are in the Argus format (extension .exp).

	Usage:
	md=setmask(md,floatingicename,groundedicename)

	Examples:
	md=setmask(md,'all','');
	md=setmask(md,'Iceshelves.exp','Islands.exp');
	"""

	elementonfloatingice = FlagElements( md, floatingicename)
	elementongroundedice = FlagElements( md, groundedicename)

	elementonfloatingice = convert( Vector{Float64}, elementonfloatingice .&  .~elementongroundedice)
	elementongroundedice = convert( Vector{Float64}, elementonfloatingice .== 0.)

	vertexonfloatingice=zeros(md.mesh.numberofvertices)
	vertexongroundedice=zeros(md.mesh.numberofvertices)

	vertexongroundedice[md.mesh.elements[findall(elementongroundedice .>0 ),:]] .= 1. # this could be an issue, since findall gives CartesianIndex, this will only take the first column from md.mesh.elements
	vertexonfloatingice[findall(vertexongroundedice .== 0.)] .= 1.

	#define levelsets
	md.mask.ocean_levelset = vertexongroundedice
	md.mask.ocean_levelset[findall(vertexongroundedice .==0.)] .= -1.
	md.mask.ice_levelset = -1*ones(md.mesh.numberofvertices)

	return md
end
#}}}
function FlagElements(md::model,region::String) #{{{

	if isempty(region)
		flags = zeros(Bool, md.mesh.numberofelements)
	elseif region == "all"
		flags = ones(Bool, md.mesh.numberofelements)
	else

		if false
			#test center of elements only
			xcenter = md.mesh.x[md.mesh.elements]*[1;1;1]/3
			ycenter = md.mesh.y[md.mesh.elements]*[1;1;1]/3
			flags = ContourToNodes(xcenter, ycenter, region, 2.)
		else
			#Follow ISSM's philosophy, all nodes need to be inside
			flags_nodes = ContourToNodes(md.mesh.x, md.mesh.y, region, 2.)
			flags = Vector{Bool}(vec((sum(flags_nodes[md.mesh.elements], dims=2).==3)))
		end
	end

	return flags
end
#}}}
function SetIceShelfBC(md::model) # {{{
	nodeonicefront=BitArray(zeros(Bool,md.mesh.numberofvertices))
	SetIceShelfBC(md, nodeonicefront)
end# }}}
function SetIceShelfBC(md::model, icefrontfile::String) # {{{
	nodeinsideicefront = ContourToNodes(md.mesh.x,md.mesh.y,icefrontfile,2.)
   nodeonicefront = (md.mesh.vertexonboundary .& nodeinsideicefront)
	SetIceShelfBC(md, nodeonicefront)
end# }}}
function SetIceShelfBC(md::model, nodeonicefront::BitVector) #{{{
	if length(nodeonicefront) == md.mesh.numberofvertices
		pos=findall(md.mesh.vertexonboundary .& .!(nodeonicefront))
		md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices)
		md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices)
		#md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices)

		# Ice front position
		md.mask.ice_levelset[findall((nodeonicefront))] .= 0.0

		if typeof(md.mesh)==DJUICE.Mesh2dTriangle
			numbernodesfront = 2
		elseif typeof(md.mesh)==DJUICE.Mesh3dPrism
			numbernodesfront = 4
		else
			error("mesh type not supported yet")
		end
		segmentsfront=md.mask.ice_levelset[md.mesh.segments[:,1:numbernodesfront]] .== 0
		segments=findall(sum(segmentsfront,dims=2).!=numbernodesfront)
		pos=md.mesh.segments[first.(Tuple.(segments)), 1:numbernodesfront]
		md.stressbalance.spcvx[pos[:]] .= 0.0
		md.stressbalance.spcvy[pos[:]] .= 0.0

		# Dirichlet Values
		if length(md.inversion.vx_obs)==md.mesh.numberofvertices & length(md.inversion.vy_obs)==md.mesh.numberofvertices
			@printf "   boundary conditions for stressbalance model: spc set as observed velocities\n"
			md.stressbalance.spcvx[pos]=md.inversion.vx_obs[pos]
			md.stressbalance.spcvy[pos]=md.inversion.vy_obs[pos]
		else
			@printf "   boundary conditions for stressbalance model: spc set as zero\n"
		end
	else
		error("Size of the input vector is not consistent with the mesh")
	end
	return md
end# }}}
