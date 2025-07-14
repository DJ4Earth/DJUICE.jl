#Contour class definition
mutable struct Contour#{{{
	id::Int64
	nods::Int64
	x::Vector{Float64}
	y::Vector{Float64}
	closed::Bool
end# }}}

