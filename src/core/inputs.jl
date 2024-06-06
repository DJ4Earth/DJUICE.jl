#Input class definition
mutable struct Input#{{{
	enum::IssmEnum
	interp::IssmEnum
	values::Vector{Float64}
	element_values::Vector{Float64}
end# }}}

#Inputs dataset definition
mutable struct Inputs #{{{
	numberofelements::Int64
	numberofvertices::Int64
	lookup::Dict{IssmEnum,Input}
end# }}}

#Inputs functions
function DuplicateInput(inputs::Inputs, old::IssmEnum, new::IssmEnum)#{{{

	#Fetch input that needs to be copied
	oldinput = inputs.lookup[old]

	if typeof(oldinput)==Input
		inputs.lookup[new] = Input(new, oldinput.interp, copy(oldinput.values), copy(oldinput.element_values))
	end

	return nothing
end#}}}
function GetInput(inputs::Inputs,enum::IssmEnum) #{{{

	#Does this input exist
	if !haskey(inputs.lookup,enum)
		error("Input ",enum," not found")
	end

	#return input
	return inputs.lookup[enum]

end#}}}
function GetInputAverageValue(input::Input) #{{{

	numnodes = NumberofNodesTria(input.interp)
	value = 0.0

	for i in 1:numnodes
		value+=input.element_values[i]
	end

	return value/numnodes

end#}}}
function GetInputDerivativeValue(input::Input,xyz_list::Matrix{Float64},gauss::GaussTria,i::Int64) #{{{

	#Get nodal function derivatives in reference element
	numnodes = NumberofNodesTria(input.interp)
	dbasis_ref = Matrix{Float64}(undef,numnodes,2)
	NodalFunctionsDerivativesReferenceTria(dbasis_ref,gauss,input.interp)

	#Get invert of the Jacobian
	Jinv = JacobianInvert(xyz_list,gauss)

	#Build dbasis:
	#[ dNi/dx ] = Jinv * [dNhat_i/dr]
	#[ dNi/dy ] =        [dNhat_i/ds]
	dbasis = Matrix{Float64}(undef,numnodes,2)
	for i in 1:3
		dbasis[i,1] = Jinv[1,1]*dbasis_ref[i,1]+Jinv[1,2]*dbasis_ref[i,2]
		dbasis[i,2] = Jinv[2,1]*dbasis_ref[i,1]+Jinv[2,2]*dbasis_ref[i,2]
	end

	#Get derivatives: dp/dx dp/dy
	dp = [0.0;0.0]
	for i in 1:3
		dp[1] += dbasis[i,1]*input.element_values[i]
		dp[2] += dbasis[i,2]*input.element_values[i]
	end

	return dp

end#}}}
function GetInputMax(input::Input) #{{{

	return maximum(input.element_values)

end#}}}
function GetInputMin(input::Input) #{{{

	return minimum(input.element_values)

end#}}}
function GetInputValue(input::Input,gauss::GaussTria,i::Int64) #{{{

	if input.interp==P0Enum
		return input.element_values[1]
	elseif input.interp==P1Enum
		value = input.element_values[1]*gauss.coords1[i] +  input.element_values[2]*gauss.coords2[i] +  input.element_values[3]*gauss.coords3[i]
	else
		error("not implemented yet")
	end

	return value

end#}}}
function SetTriaInput(inputs::Inputs,enum::IssmEnum,interp::IssmEnum,index::Int64,value::Float64) #{{{

	#Does this input exist
	if !haskey(inputs.lookup,enum)
		#it does not exist yet, we need to create it...
		@assert inputs.numberofelements > 0
		if interp==P0Enum
			inputs.lookup[enum] = Input(enum,interp,zeros(inputs.numberofelements),Vector{Float64}(undef,1))
		elseif interp==P1Enum
			inputs.lookup[enum] = Input(enum,interp,zeros(inputs.numberofvertices),Vector{Float64}(undef,3))
		else
			error("not supported yet")
		end
	end

	#Get this input and check type
	input::Input = inputs.lookup[enum]
	if typeof(input)!=Input error("input type not consistent") end
	if interp!=input.interp        error("input interpolations not consistent") end

	#set value
	input.values[index] = value

	return nothing
end#}}}
function SetTriaInput(inputs::Inputs,enum::IssmEnum,interp::IssmEnum,indices::Vector{Int64},values::Vector{Float64}) #{{{

	#Does this input exist
	if !haskey(inputs.lookup,enum)
		#it does not exist yet, we need to create it...
		@assert inputs.numberofvertices>0
		if interp==P1Enum
			inputs.lookup[enum] = Input(enum,interp,zeros(inputs.numberofvertices),Vector{Float64}(undef,3))
		else
			error("not supported yet")
		end
	end

	#Get this input and check type
	input::Input = inputs.lookup[enum]
	if typeof(input)!=Input error("input type not consistent") end
	if interp!=input.interp        error("input interpolations not consistent") end

	#set value
	input.values[indices] = values

	return nothing
end#}}}
