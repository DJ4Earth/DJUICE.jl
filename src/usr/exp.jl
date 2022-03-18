#for exptool, look into this http://juliaplots.org/MakieReferenceImages/gallery//mouse_picking/index.html

#exp object definition, constructor, and disp
mutable struct ExpStruct #{{{
	name::String
	nods::Int32
	density::Float64
	x::Vector{Float64}
   y::Vector{Float64}
	closed::Bool
end  #}}}
function ExpStruct() #{{{
	return ExpStruct("",0, 0., Vector{Float64}(undef,0), Vector{Float64}(undef,0), false)
end# }}}
function Base.show(io::IO, exp::ExpStruct)# {{{

	compact = get(io, :compact, false)

	println(io,"ExpStruct:")
	for name in fieldnames(typeof(exp))
		a=getfield(exp,name)
		print(io,"   $(name) = ")
		if !isempty(a)
			if compact && eltype(a)<:Number && length(a)>3
				println(io, typeof(a), " of size ", size(a))
			else
				println(io,a)
			end
		else
			println(io,"empty")
		end
	end
end# }}}

#methods
#expread {{{
"""
	EXPREAD - read a file exp and build a Structure

	This function reads an *.exp* and builds a structure containing the fields x
	and y corresponding to the coordinates, one for the filename of the exp
	file, for the density, for the nodes, and a field 'closed' to indicate if the
	domain is closed.

	Usage:
		exp=expread(filename)

# Examples:
```julia-repl
julia> exp=expread('domainoutline.exp')
```

# Arguments:
- filename: the ARGUS file to read
"""
function expread(filename::String)

	#initialize output
	contours = Vector{ExpStruct}(undef, 0)

	#open file
	f = open(filename, "r") do f

		#initialize some variables
		nprof = 0
		line = 1

		while !eof(f)

			#read first line
			A = readline(f); line += 1

			#if isempty, go to the next line and try again
			if isempty(A)
				continue
			else
				#initialize new profile
				nprof += 1 
				exp = ExpStruct();
			end

			#extract profile name
			if A[1:8]!="## Name:"
				println("line $(line): $(A)")
				error("Unexpected exp file formatting") 
			end
			if length(A)>8
				exp.name = A[9:end]
			end

			#read Icon
			A = readline(f); line += 1
			if A[1:8]!="## Icon:" error("Unexpected exp file formatting") end

			#read Info
			A = readline(f); line += 1
			if A[1:14]!="# Points Count"
				println("line $(line): $(A)")
				error("Unexpected exp file formatting") 
			end

			#Reads number of nods and density
			A = readline(f); line += 1
			A = parse.(Float64, split(A))
			if length(A) != 2 error("Unexpected exp file formatting") end
			exp.nods = A[1]; exp.density = A[2]

			#Allocate arrays
			if exp.nods<=0 error("Unexpected exp file formatting") end
			exp.x = Vector{Float64}(undef,exp.nods)
			exp.y = Vector{Float64}(undef,exp.nods)

			#Read coordinates
			A = readline(f); line += 1
			if A[1:13]!="# X pos Y pos" error("Unexpected exp file formatting") end
			for i in 1:exp.nods
				A = readline(f); line += 1
				A = parse.(Float64, split(A))
				if length(A) != 2 error("Unexpected exp file formatting") end
				if any(isnan.(A)) error("NaNs found in coordinate") end
				exp.x[i] = A[1]; exp.y[i] = A[2]
			end

			#check if closed
			if exp.nods>1 && exp.x[1]==exp.x[end] && exp.y[1]==exp.y[end]
				exp.closed = true
			else
				exp.closed = false
			end

			#add profile to list
			push!(contours, exp)
		end
	end

	return contours
end# }}}
#ContourToNodes{{{
"""
	ContourToNodes - Flag points that are in contour

   More doc to come later....

	Usage:
		exp=expread(filename)

# Examples:
```julia-repl
julia> exp=expread('domainoutline.exp')
```

# Arguments:
- filename: the ARGUS file to read
"""
function ContourToNodes(x::Vector{Float64},y::Vector{Float64},filename::String,edgevalue::Float64)

	#Read input file
	contours = expread(filename)

	#Initialize output
	nbpts = length(x)
	flags = zeros(Bool,nbpts)

	#Loop over contours
	for c in 1:length(contours)

		#Get current contours
		contour = contours[c]
		xp      = contour.x
		yp      = contour.y

		#Check that we are within box
		xmin = minimum(xp); xmax = maximum(xp)
		ymin = minimum(yp); ymax = maximum(yp)

		#Loop over all points provided
		for ii in 1:nbpts

			#If this node is already within one of the contours, do not change it
			if(flags[ii]) continue end

			#Are we within bounds?
			if(x[ii]<xmin || x[ii]>xmax || y[ii]<ymin || y[ii]>ymax) continue end

			#we are potentially inside... perform pnpoly test
			flags[ii] = pnpoly(xp, yp, x[ii], y[ii], edgevalue)
		end
	end
	
	return flags
end# }}}

function pnpoly(xp::Vector{Float64},yp::Vector{Float64},x::Float64,y::Float64,edgevalue::Float64) #{{{

	npol = length(xp)

	#Do we need to test for colinearity?
	if(edgevalue!=2)
		i = 1
		j = npol
		while(i<=npol)

			n1 = (yp[i]-yp[j])^2 + (xp[i]-xp[j])^2
			n2 = (y-yp[j])^2 + (x-xp[j])^2

			normp=sqrt(n1*n2)
			scalar=(yp[i]-yp[j])*(y-yp[j])+(xp[i]-xp[j])*(x-xp[j])

			if (scalar == normp)
				if (n2<=n1)
					return edgevalue
				end
			end

			j =  i
			i += 1
		end
	end

	#second test : point is neither on a vertex, nor on a side, where is it ?
	i = 1
	j = npol
	c = false
	while(i<=npol)
		if (((yp[i]<=y && y<yp[j]) || (yp[j]<=y && y<yp[i])) &&
            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
         c = !c
		end

		j =  i
		i += 1
	end
	return c
end# }}}
