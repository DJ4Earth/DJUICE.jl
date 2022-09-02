#utils
function issmdir() #{{{
	issmdir = ENV["ISSM_DIR"]

	if isempty(issmdir)
		error("Could not determine the location of ISSM")
	else
		return issmdir
	end
end#}}}
function archread(filename::String,variablename::String) #{{{

	#initialize variables
	found = false

	#open file
	output = open(filename, "r") do f

		while !eof(f)
			reclen  = bswap(read(f, Int32))
			rectype = bswap(read(f, Int32))
			if rectype!=1
				error("Expected variable of type string")
			else
				fieldname_length = bswap(read(f, Int32))
				field_name = String(read(f, fieldname_length))
			end
			rec_length = bswap(read(f, Int32))
			field_type = bswap(read(f, Int32))
			if field_type==2
				data = bswap(read(f, Float64))
			elseif field_type==3
				rows = bswap(read(f, Int32))
				cols = bswap(read(f, Int32))
				data = reinterpret(Float64, read(f, sizeof(Float64)*rows*cols))
				data .= ntoh.(data)
				data = reshape(data, (rows,cols))
				data = collect(data)
				if cols == 1
					data = vec(data)
				end
			else
				error("Error: Encountered invalid field type when reading data.")
			end

			if field_name == variablename
				found = true
				return data
			end
		end
	end

	return output
end# }}}
function InterpFromMeshToMesh2d(index_data::Array,x_data::Vector,y_data::Vector,data::Vector,xout::Vector,yout::Vector,default::Float64=NaN) #{{{

	#Allocate output
	nods_out = length(xout)
	data_out = default*ones(nods_out)

	#Interpolation type
	data_length = size(data,1)
	nods_data   = length(x_data)
	nels_data   = size(index_data,1)
	if(data_length==nods_data)
		interpolation_type=1;
	elseif (data_length==nels_data)
		interpolation_type=2
	else
		error("length of vector data not supported yet. It should be of length (number of nodes) or (number of elements)!")
	end
	xmin = minimum(xout); xmax = maximum(xout)
	ymin = minimum(yout); ymax = maximum(yout)

	for i in 1:nels_data

		#skip element if no overlap
		if (minimum(x_data[index_data[i,:]]) > xmax) continue end
		if (minimum(y_data[index_data[i,:]]) > ymax) continue end
		if (maximum(x_data[index_data[i,:]]) < xmin) continue end
		if (maximum(y_data[index_data[i,:]]) < ymin) continue end

		#get area of the current element (Jacobian = 2 * area)*/
		#area =x2 * y3 - y2*x3 + x1 * y2 - y1 * x2 + x3 * y1 - y3 * x1;
		area = (x_data[index_data[i,2]]*y_data[index_data[i,3]]-y_data[index_data[i,2]]*x_data[index_data[i,3]] 
				  +  x_data[index_data[i,1]]*y_data[index_data[i,2]]-y_data[index_data[i,1]]*x_data[index_data[i,2]] 
				  +  x_data[index_data[i,3]]*y_data[index_data[i,1]]-y_data[index_data[i,3]]*x_data[index_data[i,1]])

		for j in 1:nods_out
			#Get first area coordinate = det(x-x3  x2-x3 ; y-y3   y2-y3)/area
			area_1=((xout[j]-x_data[index_data[i,3]])*(y_data[index_data[i,2]]-y_data[index_data[i,3]])
					 -  (yout[j]-y_data[index_data[i,3]])*(x_data[index_data[i,2]]-x_data[index_data[i,3]]))/area
			#Get second area coordinate =det(x1-x3  x-x3 ; y1-y3   y-y3)/area
			area_2=((x_data[index_data[i,1]]-x_data[index_data[i,3]])*(yout[j]-y_data[index_data[i,3]])
					  - (y_data[index_data[i,1]]-y_data[index_data[i,3]])*(xout[j]-x_data[index_data[i,3]]))/area
			#Get third area coordinate = 1-area1-area2
			area_3=1-area_1-area_2

			if (area_1>=0 && area_2>=0 && area_3>=0)
				if (interpolation_type==1)
					#nodal interpolation
					data_out[j]=area_1*data[index_data[i,1]]+area_2*data[index_data[i,2]]+area_3*data[index_data[i,3]];
				else
					#element interpolation
					data_out[j]=data[i];
				end
			end
		end
	end
	return data_out

	#OLD STUFF!!! not working...

	#prepare input arrays
	nods = Cint(length(x))
	nels = Cint(size(index,1))
	nods_interp = Cint(length(xout))
	Cindex=Array{Cint,1}(undef,length(index))
	for i in 1:size(index,1)
		for j in 1:3
			Cindex[(i-1)*3+j] = Int32(index[i,j])
		end
	end
	Cx    = Array{Cdouble,1}(undef,nods)
	Cy    = Array{Cdouble,1}(undef,nods)
	Cdata = Array{Cdouble,1}(undef,nods)
	for i in 1:nods
		Cx[i]    = x[i]
		Cy[i]    = y[i]
		Cdata[i] = data[i]
	end
	Cxout = Array{Cdouble,1}(undef,nods_interp)
	Cyout = Array{Cdouble,1}(undef,nods_interp)
	for i in 1:nods_interp
		Cxout[i] = xout[i]
		Cyout[i] = yout[i]
	end

	Cdataout = Vector{Float64}(undef,nods_interp)

	#This is not working....
	rc=ccall( (:InterpFromMeshToMesh2dx,"libISSMCore"),
				Cint, (Ptr{Ptr{Cdouble}},Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint),
				Ref(Ref(Cdataout)), Ref(Cindex), Ref(Cx), Ref(Cy), nods, nels,
				Ref(Cdata), nods, 1, Ref(Cxout), Ref(Cyout), nods_interp)

	#Process output
	dataout = Vector{Float64}(undef,nods_interp)
	for i in 1:nods_interp
		dataout[i] = Cdataout[i]
	end

	return dataout
end #}}}
function InterpFromMeshToMesh2d2(index_data::Array,x_data::Vector,y_data::Vector,data::Vector,xout::Vector,yout::Vector) #{{{

	#prepare input arrays
	nods = Cint(length(x_data))
	nels = Cint(size(index_data,1))
	nods_interp = Cint(length(xout))
	Cindex=Array{Cint,1}(undef,length(index_data))
	for i in 1:size(index_data,1)
		for j in 1:3
			Cindex[(i-1)*3+j] = Int32(index_data[i,j])
		end
	end
	Cx    = Array{Cdouble,1}(undef,nods)
	Cy    = Array{Cdouble,1}(undef,nods)
	Cdata = Array{Cdouble,1}(undef,nods)
	for i in 1:nods
		Cx[i]    = x_data[i]
		Cy[i]    = y_data[i]
		Cdata[i] = data[i]
	end
	Cxout    = Array{Cdouble,1}(undef,nods_interp)
	Cyout    = Array{Cdouble,1}(undef,nods_interp)
	Cdataout = Array{Cdouble,1}(undef,nods_interp)
	for i in 1:nods_interp
		Cxout[i] = xout[i]
		Cyout[i] = yout[i]
	end

	#This is not working....
	#rc=ccall( (:InterpFromMeshToMesh2dx,"../bamg/libBamg.so"),
	#			Cint, (Ptr{Cdouble},Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint),
	#			Ref(Cdataout), Ref(Cindex), Ref(Cx), Ref(Cy), nods, nels,
	#			Ref(Cdata), nods, 1, Ref(Cxout), Ref(Cyout), nods_interp)
	#rc=ccall( (:InterpFromMeshToMesh2dx,"../bamg/libBamg.so"),
	#			Cint, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint),
	#			Ref(Cindex), Ref(Cx), Ref(Cy), nods, nels)
	#
	#
	dataout = Vector{Float64}(undef,nods_interp)
	rc=ccall( (:InterpFromMeshToMesh2dx3,"/Users/mmorligh/Desktop/issmuci/trunk-jpl/src/jl/bamg/libBamg.dylib"),
				Cint, (Ptr{Cdouble}, Cint),
				dataout, nods_interp)

	#Process output
	for i in 1:nods_interp
		dataout[i] = Cdataout[i]
	end

	return dataout
end #}}}
function IssmStructDisp(io::IO, modelfield::Any) # {{{
	println(io,typeof(modelfield),":")
	for name in fieldnames(typeof(modelfield))
		a=getfield(modelfield,name)
		#print(io,"   $(name) = ")
		@printf "%19s: " name
		if isa(a,String)
			println(io, a)
		elseif isa(a, Flux.Chain)
			println(io, "Flux.", a)
		elseif length(a)>1
			if !isempty(a)
				println(io, typeof(a), " of size ", size(a))
			else
				println(io,"empty")
			end
		else
			println(io, a)
		end
	end
end #}}}
function meshgrid(x::Vector, y::Vector) # {{{
	X = [i for i in x, j in 1:length(y)]
	Y = [j for i in 1:length(x), j in y]
	return X, Y
end #}}}
