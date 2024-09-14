#Input class definition
mutable struct Input#{{{
	enum::IssmEnum
	interp::IssmEnum
	values::Vector{Float64}
	element_values::Vector{Float64}
end# }}}
mutable struct TransientInput #{{{
	enum::IssmEnum
	times::Vector{Float64}
	N::Int64
	inputs::Vector{Input}
	parameters::Parameters
	currentinput::Input
	currentstep::Float64
end# }}}

#Inputs dataset definition
mutable struct Inputs #{{{
	numberofelements::Int64
	numberofvertices::Int64
	lookup::Dict{IssmEnum,Union{Input,TransientInput}}
end# }}}

#Inputs functions
function Configure(inputs::Inputs, parameters::Parameters) #{{{

	#Loop over all inputs and find transient inputs
	for (key, value) in inputs.lookup
		if isa(value, TransientInput)
			value.parameters = parameters
		end
	end

	return nothing

end#}}}
function DuplicateInput(inputs::Inputs, old::IssmEnum, new::IssmEnum)#{{{

	#Fetch input that needs to be copied
	oldinput = inputs.lookup[old]

	if typeof(oldinput)==Input
		inputs.lookup[new] = Input(new, oldinput.interp, copy(oldinput.values), copy(oldinput.element_values))
	end

	return nothing
end#}}}
function GetInput(inputs::Inputs, enum::IssmEnum) #{{{

	#Does this input exist
	if !haskey(inputs.lookup,enum)
		error("Input ",enum," not found")
	end

	#return input
	input = inputs.lookup[enum]
	if isa(input, TransientInput)

		#Find time in parameters
		time = FindParam(Float64, input.parameters, TimeEnum)

		#Interpolate input if necessary
		SetCurrentTimeInput(input, time)
		input = input.currentinput
	end

	#Make sure this is a regular input
	if ~isa(input, Input) error("Input ",enum," is not an Input") end
	return input

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
	if interp!=input.interp error("input interpolations not consistent") end

	#set value
	input.values[indices] = values

	return nothing
end#}}}

#TransientInput functions
function GetTransientInput(inputs::Inputs, enum::IssmEnum) #{{{

	#Does this input exist
	if !haskey(inputs.lookup,enum)
		error("Input ",enum," not found")
	end

	#return input
	if typeof(inputs.lookup[enum])!=TransientInput error("Input ",enum," is not a TransientInput") end
	return inputs.lookup[enum]

end#}}}
function AddTimeInput(inputs::Inputs, tr_input::TransientInput, index::Int64, interp::IssmEnum,indices::Vector{Int64},values::Vector{Float64}) #{{{

	#Check index
	@assert index>0 && index<=tr_input.N

	#Is input already assigned?
	if ~isassigned(tr_input.inputs, index)
		@assert inputs.numberofelements > 0
		if interp==P1Enum
			tr_input.inputs[index] = Input(tr_input.enum, interp,zeros(inputs.numberofvertices), Vector{Float64}(undef,3))
		else
			error("not supported yet")
		end
	end

	#set value
	tr_input.inputs[index].values[indices] = values

	return nothing
end#}}}
function SetTransientInput(inputs::Inputs,enum::IssmEnum,times::Vector{Float64}) #{{{

	#Does this input exist?
	if !haskey(inputs.lookup, enum)
		#it does not exist yet, we need to create it...
		N = length(times)
		@assert N>0
		inputs.lookup[enum] = TransientInput(enum, times, length(times), Vector{Input}(undef,N), Parameters(), Input(enum, P0Enum ,zeros(0), Vector{Float64}(undef,0)), -1.)
	end

	#Some checks that everything is consistent
	transientinput::TransientInput = inputs.lookup[enum]
	if typeof(transientinput)!=TransientInput error("input type not consistent") end
	if length(times)!=transientinput.N        error("input time series consistent") end

	return nothing
end#}}}
function SetCurrentTimeInput(tr_input::TransientInput, time::Float64) #{{{

	#Binary search where time is in our series
	p = searchsorted(tr_input.times, time)

	println("time is ",time, " and p is ", p)

	if p.stop==0
		#Before our time series
		@assert time<=tr_input.times[1]

		#If already processed return
		if(tr_input.currentstep==0.) return nothing end

		#Copy first input
		tr_input.currentinput = tr_input.inputs[0]
		tr_input.currentstep  = 0.

	elseif p.start==length(tr_input.times)+1
		#After the end of our time series
		@assert time>=tr_input.times[end]

		#If already processed return
		if(tr_input.currentstep==Float64(tr_input.N)) return nothing end

		#Copy last input
		tr_input.currentinput = tr_input.inputs[end]
		tr_input.currentstep  = Float64(tr_input.N)
	else
		#General case...
		error("Not implemented yet, see C++ code")
	end

	return nothing
end#}}}
