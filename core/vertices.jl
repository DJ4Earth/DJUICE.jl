#Vertex class definition
mutable struct Vertex#{{{
	sid::Int64
	x::Float64
	y::Float64
	z::Float64
end# }}}

#vertices functions
function GetVerticesCoordinates(vertices::Vector{Vertex}) #{{{

	#Intermediaries
	nbv = length(vertices)

	#Allocate
	xyz_list = Matrix{Float64}(undef,nbv,3)

	#Assign value to xyz_list
	for i in 1:nbv
		xyz_list[i,1]=vertices[i].x
		xyz_list[i,2]=vertices[i].y
		xyz_list[i,3]=vertices[i].z
	end

	return xyz_list
end #}}}
