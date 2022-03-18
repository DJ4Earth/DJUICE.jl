#Constraint class definition
struct Constraint #{{{
	id::Int64
	nodeid::Int64
	dof::Int8
	value::Float64
end# }}}

#Constraint functions
function ConstrainNode(constraint::Constraint,nodes::Vector{Node},parameters::Parameters) #{{{

	#Chase through nodes and find the node to which this SpcStatic apply
	node = nodes[constraint.nodeid]

	#Apply Constraint
	ApplyConstraint(node, constraint.dof, constraint.value)

	return nothing
end# }}}
