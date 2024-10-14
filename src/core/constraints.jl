#Constraint class definition
abstract type AbstractConstraint end
struct Constraint <: AbstractConstraint #{{{
	id::Int64
	nodeid::Int64
	dof::Int8
	value::Float64
end# }}}
struct ConstraintTransient <: AbstractConstraint #{{{
	id::Int64
	nodeid::Int64
	dof::Int8
	times::Vector{Float64}
	values::Vector{Float64}
end# }}}

#Constraint functions
function ConstrainNode(constraint::Constraint,nodes::Vector{Node},parameters::Parameters) #{{{

	#Chase through nodes and find the node to which this SpcStatic applies
	node = nodes[constraint.nodeid]

	#Apply Constraint
	ApplyConstraint(node, constraint.dof, constraint.value)

	return nothing
end# }}}
function ConstrainNode(constraint::ConstraintTransient, nodes::Vector{Node}, parameters::Parameters) #{{{

	#Chase through nodes and find the node to which this SpcTransient applies
	node = nodes[constraint.nodeid]

	#Find time in parameters
	time = FindParam(Float64, parameters, TimeEnum)

	# Now, go fetch value for this time:
	if (time<=times[1])
		value = values[1]
		found = true
	elseif (time>=times[end])
		value = values[end]
		found = true
	else
		for i in 1:length(times)-1
			if (times[i]<=time && time<times[i+1])
				alpha = (time-times[i])/(times[i+1]-times[i]);
				@assert alpha>=0.0 && alpha<=1.0
				value = (1-alpha)*values[i] + alpha*values[i+1];
				found = true
				break
			end
		end
	end

	if(!found) error("could not find time segment for constraint") end

	#Apply Constraint
	if isnan(value)
		RelaxConstraint(node, constraint.dof)
	else
		ApplyConstraint(node, constraint.dof, constraint.value)
	end

	return nothing
end# }}}
