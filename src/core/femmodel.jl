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
	segments = Vector{Float64}(undef,0)

	for i in 1:length(femmodel.elements)
		WriteFieldIsovalueSegment!(femmodel.elements[i], segments, fieldnum, fieldvalue)
	end

	# TODO: parallel version
	
	# Add distance input to all elements

	return nothing
end#}}}
