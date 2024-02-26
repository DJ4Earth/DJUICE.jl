using Triangulate

"""
TRIANGLE - create model mesh using the triangle package

	This function creates a model mesh using Triangle and a domain outline, to
	within a certain resolution
#Arguments
 - md is a model tuple
 - domainname is the name of an Argus domain outline file
 - resolution:  is a characteristic length for the mesh (same unit as the domain outline unit)

# Usage:
 - md=triangle(md,domainname,resolution)
# Examples:
 - md=triangle(md,'DomainOutline.exp',1000);
 - md=triangle(md,'DomainOutline.exp','Rifts.exp',1500);
"""
function triangle(md::model,domainname::String,resolution::Float64) #{{{

	#read input file
	contours = expread(domainname)
	return triangle(md, contours, resolution)
end#}}}

function triangle(md::model,contours::Vector{ExpStruct},resolution::Float64) #{{{
	area     = resolution^2

	#Initialize input structure
	triin=Triangulate.TriangulateIO()

	#Construct input structure
	numberofpoints   = 0
	numberofsegments = 0
	for i in 1:length(contours)
		numberofpoints   += contours[i].nods-1
		numberofsegments += contours[i].nods-1
	end
	numberofpointattributes = 1

	pointlist=Array{Cdouble,2}(undef,2,numberofpoints)
	count = 0
	for i in 1:length(contours)
		nods = contours[i].nods
		pointlist[1,count+1:count+nods-1] = contours[i].x[1:end-1]
		pointlist[2,count+1:count+nods-1] = contours[i].y[1:end-1]
		count += (nods-1)
	end
	pointattributelist=Array{Cdouble,1}(undef,numberofpoints)
	pointmarkerlist=Array{Cint,1}(undef,numberofpoints)
	for i in 1:numberofpoints
		pointmarkerlist[i]=0
		pointattributelist[i]=0.
	end

	counter=0;
   backcounter=0;
	segmentlist=Array{Cint,2}(undef,2,numberofsegments)
	segmentmarkerlist=Array{Cint,1}(undef,numberofsegments)
	segmentmarkerlist[:].=0
	for i in 1:length(contours)
		nods = contours[i].nods
		segmentlist[1,counter+1:counter+nods-2] = collect(counter+0:counter+nods-3)
		segmentlist[2,counter+1:counter+nods-2] = collect(counter+1:counter+nods-2)
		counter+=nods-2
		#close profile
		segmentlist[1,counter+1]=counter
		segmentlist[2,counter+1]=backcounter
		counter+=1
		backcounter=counter
	end

	numberofregions = 0
	numberofholes = length(contours)-1
	holelist = Array{Cdouble,2}(undef,2,numberofholes)
	if numberofholes>0
		 for i in 2:length(contours)
			 xA=contours[i].x[1]; xB=contours[i].x[end-1]
			 yA=contours[i].y[1]; yB=contours[i].y[end-1]
			 xC=(xA+xB)/2;        yC=(yA+yB)/2;
			 xD=xC+tan(10. /180. *pi)*(yC-yA);
			 yD=yC+tan(10. /180. *pi)*(xA-xC);
			 xE=xC-tan(10. /180. *pi)*(yC-yA);
			 yE=yC-tan(10. /180. *pi)*(xA-xC);
			 holelist[1,i-1] = xD
			 holelist[2,i-1] = yD
		 end
	end

	#based on this, prepare input structure
	triin.pointlist=pointlist
	triin.segmentlist=segmentlist
	triin.segmentmarkerlist = segmentmarkerlist
	triin.holelist=holelist

	#Triangulate
	triangle_switches = "pQzDq30ia"*@sprintf("%lf",area) #replace V by Q to quiet down the logging
	(triout, vorout)=triangulate(triangle_switches, triin)

	#assign output
	md.mesh = Mesh2dTriangle()
	md.mesh.numberofvertices = size(triout.pointlist, 2)
	md.mesh.numberofelements = size(triout.trianglelist, 2)
	md.mesh.x                = triout.pointlist[1,:]
	md.mesh.y                = triout.pointlist[2,:]
	md.mesh.elements         = convert(Matrix{Int64}, triout.trianglelist' .+ 1)
	md.mesh.segments         = collect(triout.segmentlist' .+ 1)

	#post processing
	md.mesh.vertexonboundary = zeros(Bool,md.mesh.numberofvertices)
	md.mesh.vertexonboundary[md.mesh.segments] .= true

	return md
end#}}}
