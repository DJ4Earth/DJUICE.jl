#Parameter class definition
abstract type Parameter end
struct DoubleParam <: Parameter #{{{
	enum::IssmEnum
	value::Float64
end# }}}
struct IntParam <: Parameter #{{{
	enum::IssmEnum
	value::Int64
end# }}}
struct BoolParam <: Parameter #{{{
	enum::IssmEnum
	value::Bool
end# }}}

#Parameters dataset class definition
mutable struct Parameters #{{{
	lookup::Dict{IssmEnum,Parameter}
	double_lookup::Dict{IssmEnum,DoubleParam}
end# }}}

#Parameter functions
function GetParameterValue(param::DoubleParam) #{{{
	return param.value
end#}}}
function GetParameterValue(param::IntParam) #{{{
	return param.value
end#}}}
function GetParameterValue(param::BoolParam) #{{{
	return param.value
end#}}}

#Parameters functions
function AddParam(parameters::Parameters,value::Float64,enum::IssmEnum) #{{{

	parameters.double_lookup[enum] = DoubleParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Int64, enum::IssmEnum) #{{{

	parameters.lookup[enum] = IntParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Bool, enum::IssmEnum) #{{{

	parameters.lookup[enum] = BoolParam(enum,value)

	return nothing
end#}}}

@noinline function FindParam(::Type{Float64}, parameters::Parameters,enum::IssmEnum) #{{{

	param = parameters.double_lookup[enum]
	return GetParameterValue(param)::Float64

end#}}}

@noinline function FindParam(::Type{T}, parameters::Parameters,enum::IssmEnum) where T #{{{

	param = parameters.lookup[enum]
	return GetParameterValue(param)::T

end#}}}
