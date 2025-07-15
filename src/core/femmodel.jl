#femmodel class definition
mutable struct FemModel #{{{
	analyses::Vector{Analysis}
	elements::Vector{Tria}
	vertices::Vector{Vertex}

	nodes::Vector{Node}
	nodes_list::Vector{Vector{Node}}

	parameters::Parameters
	inputs::Inputs

	constraints::Vector{AbstractConstraint}
	constraints_list::Vector{Vector{AbstractConstraint}}

	#loads::Vector{Loads}

	results::Vector{Result}
end#}}}

#femmodel functions
function SetCurrentConfiguration!(femmodel::FemModel, analysis::Analysis) #{{{

	#Find the index of this analysis
	index = -1
	for i in 1:length(femmodel.analyses)
		if(typeof(femmodel.analyses[i]) == typeof(analysis)) index = i end
	end
	if(index<1) error("Could not find analysis ",analysis, " in femmodel") end

	#Plug right nodes onto element
	for i in 1:length(femmodel.elements)
		femmodel.elements[i].nodes = femmodel.elements[i].nodes_list[index]
	end

	#Plug in nodes and other datasets
	femmodel.nodes       = femmodel.nodes_list[index]
	femmodel.constraints = femmodel.constraints_list[index]

	return nothing
end#}}}
function AddResult!(femmodel::FemModel, result::Result) #{{{
	# TODO: maybe need to check if the result already exist, then overwrite
	push!(femmodel.results, result)
end#}}}
function DistanceToFieldValue!(femmodel::FemModel, fieldnum::IssmEnum, fieldvalue::Float64, distanceenum::IssmEnum) #{{{

	# get all segments
	segments = Vector{Contour}(undef,0)

	for i in 1:length(femmodel.elements)
		WriteFieldIsovalueSegment!(femmodel.elements[i], segments, fieldnum, fieldvalue)
	end

	# TODO: parallel version

	# Add distance input to all elements
	distances = Vector{Float64}(undef,length(femmodel.vertices))
	last = 0
	for i in 1:length(femmodel.vertices)
		x = femmodel.vertices[i].x
		y = femmodel.vertices[i].y
		if (last > 1)
			dmin = (segments[last].x[1]-x)^2 + (segments[last].y[1]-y)^2
		else
			dmin = 1.e50
		end

		for i in 1:length(segments)
			# 	Skip if tip is more than 10xdmin away
			if ((segments[i].x[1]-x)^2.0 > 10*dmin) || ((segments[i].y[1]-y)^2.0 > 10*dmin)
				continue
			end

			l2 = (segments[i].x[2]-segments[i].x[1])^2 + (segments[i].y[2]-segments[i].y[1])^2
			if (l2 == 0.)
				d = (x-segments[i].x[1])^2 + (y-segments[i].y[1])^2
				if (d<dmin)
					dmin = d
					last = i
				end
				continue
			end
			# Consider the line extending the segment, parameterized as v + t (w - v).
			# We find projection of point p onto the line.
			# It falls where t = [(p-v) . (w-v)] / |w-v|^2
			t = ((x-segments[i].x[1])*(segments[i].x[2]-segments[i].x[1]) + (y-segments[i].y[1])*(segments[i].y[2]-segments[i].y[1])) / l2

			if (t < 0.0)
				# Beyond the 'v' end of the segment
				d = (x-segments[i].x[1])^2 + (y-segments[i].y[1])^2
			elseif (t > 1.0)
				# Beyond the 'w' end of the segment
				d = (x-segments[i].x[2])^2 + (y-segments[i].y[2])^2
			else
				# Projection falls on the segment
				xn = segments[i].x[1] + t * (segments[i].x[2]-segments[i].x[1])
				yn = segments[i].y[1] + t * (segments[i].y[2]-segments[i].y[1])
				d = (x-xn)^2 + (y-yn)^2
			end

			if (d<dmin)
				dmin = d
				last = i
			end
		end
		distances[i] = sqrt(dmin)
	end

	# update signed distance
	for i in 1:length(femmodel.elements)
		CreateDistanceInputFromSegmentlist(femmodel.elements[i], distances, distanceenum)
	end

	return nothing
end#}}}
