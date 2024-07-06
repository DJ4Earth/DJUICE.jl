# setmask {{{
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
function setmask(md::model,floatingicename::String,groundedicename::String)

	elementonfloatingice = FlagElements( md, floatingicename)
	elementongroundedice = FlagElements( md, groundedicename)

	elementonfloatingice = convert( Vector{Float64}, elementonfloatingice .&  .~elementongroundedice)
	elementongroundedice = convert( Vector{Float64}, elementonfloatingice .== 0.)

	vertexonfloatingice=zeros(md.mesh.numberofvertices)
	vertexongroundedice=zeros(md.mesh.numberofvertices)

	vertexongroundedice[md.mesh.elements[findall(elementongroundedice .>0 ),:]] .= 1.
	vertexonfloatingice[findall(vertexongroundedice .== 0.)] .= 1.

	#define levelsets
	md.mask.ocean_levelset = vertexongroundedice
	md.mask.ocean_levelset[findall(vertexongroundedice .==0.)] .= -1.
	md.mask.ice_levelset = -1*ones(md.mesh.numberofvertices)

	return md
end
#}}}
# FlagElements {{{
function FlagElements(md::model,region::String)

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
