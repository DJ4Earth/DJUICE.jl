
#Class Triangle's triangulateio
mutable struct CTriangulateIO #{{{

    pointlist :: Ptr{Cdouble}
    pointattributelist :: Ptr{Cdouble}
    pointmarkerlist :: Ptr{Cint}
    numberofpoints :: Cint
    numberofpointattributes :: Cint

    trianglelist :: Ptr{Cint}
    triangleattributelist :: Ptr{Cdouble}
    trianglearealist :: Ptr{Cdouble}
    neighborlist :: Ptr{Cint}
    numberoftriangles :: Cint
    numberofcorners :: Cint
    numberoftriangleattributes :: Cint

    segmentlist :: Ptr{Cint}
    segmentmarkerlist :: Ptr{Cint}
    numberofsegments :: Cint

    holelist :: Ptr{Cdouble}
    numberofholes :: Cint

    regionlist :: Ptr{Cdouble}
    numberofregions :: Cint

    edgelist :: Ptr{Cint}
    edgemarkerlist :: Ptr{Cint}
    normlist :: Ptr{Cdouble}
    numberofedges :: Cint
 end  #}}}
function CTriangulateIO() #{{{
	return CTriangulateIO(C_NULL, C_NULL, C_NULL, 0, 0,
								 C_NULL, C_NULL, C_NULL, C_NULL, 0, 0, 0,
								 C_NULL, C_NULL, 0,
								 C_NULL, 0,
								 C_NULL, 0,
								 C_NULL, C_NULL, C_NULL, 0)
end# }}}
function Base.show(io::IO, tio::CTriangulateIO)# {{{
	println(io,"CTriangulateIO(")
	for name in fieldnames(typeof(tio))
		a=getfield(tio,name)
		print(io,"$(name) = ")
		println(io,a)
	end
	println(io,")")
end# }}}

using Printf #needed for sprintf

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
function triangle_issm(md::model,domainname::String,resolution::Float64) #{{{

	#read input file
	contours = expread(domainname)
	area     = resolution^2

	#Initialize i/o structures
	ctio_in  = CTriangulateIO();
	ctio_out = CTriangulateIO();
	vor_out  = CTriangulateIO();

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
	ctio_in.numberofpoints = numberofpoints
	ctio_in.pointlist=pointer(pointlist)
	ctio_in.numberofpointattributes=numberofpointattributes
	ctio_in.pointattributelist=pointer(pointattributelist)
	ctio_in.pointmarkerlist=pointer(pointmarkerlist)
	ctio_in.numberofsegments=numberofsegments
	ctio_in.segmentlist=pointer(segmentlist)
	ctio_in.segmentmarkerlist = pointer(segmentmarkerlist)
	ctio_in.numberofholes=numberofholes
	ctio_in.holelist=pointer(holelist)
	ctio_in.numberofregions=0

	#Call triangle using ISSM's default options
	triangle_switches = "pQzDq30ia"*@sprintf("%lf",area) #replace V by Q to quiet down the logging
	#rc=ccall( (:triangulate,"libtriangle"),

	try rc=ccall( (:triangulate,issmdir()*"/externalpackages/triangle/src/libtriangle."*(@static Sys.islinux() ? :"so" : (@static Sys.isapple() ? :"dylib" : :"so"))),
				Cint, ( Cstring, Ref{CTriangulateIO}, Ref{CTriangulateIO}, Ref{CTriangulateIO}),
				triangle_switches, Ref(ctio_in), Ref(ctio_out), Ref(vor_out))
	catch LoadError
		rc=ccall( (:triangulate,issmdir()*"/externalpackages/triangle/install/lib/libtriangle."*(@static Sys.islinux() ? :"so" : (@static Sys.isapple() ? :"dylib" : :"so"))),
					Cint, ( Cstring, Ref{CTriangulateIO}, Ref{CTriangulateIO}, Ref{CTriangulateIO}),
					triangle_switches, Ref(ctio_in), Ref(ctio_out), Ref(vor_out))
	end

	#post process output
	points    = convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio_out.pointlist,    (2,Int(ctio_out.numberofpoints)), own=true))'
	triangles = convert(Array{Cint,2},    Base.unsafe_wrap(Array, ctio_out.trianglelist, (3,Int(ctio_out.numberoftriangles)), own=true))' .+1
	segments  = convert(Array{Cint,2},    Base.unsafe_wrap(Array, ctio_out.segmentlist,  (2,Int(ctio_out.numberofsegments)), own=true))' .+1

	#assign output
	md.mesh = Mesh2dTriangle()
	md.mesh.numberofvertices = ctio_out.numberofpoints
	md.mesh.numberofelements = ctio_out.numberoftriangles
	md.mesh.x                = points[:,1]
	md.mesh.y                = points[:,2]
	md.mesh.elements         = triangles
	md.mesh.segments         = segments

	#post processing
	md.mesh.vertexonboundary = zeros(Bool,md.mesh.numberofvertices)
	md.mesh.vertexonboundary[md.mesh.segments] .= true

	return md
end#}}}
