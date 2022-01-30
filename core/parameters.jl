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

	parameters.lookup[enum] = DoubleParam(enum,value)

end#}}}
function AddParam(parameters::Parameters,value::Int64, enum::IssmEnum) #{{{

	parameters.lookup[enum] = IntParam(enum,value)

end#}}}
function AddParam(parameters::Parameters,value::Bool, enum::IssmEnum) #{{{

	parameters.lookup[enum] = BoolParam(enum,value)

end#}}}
function FindParam(::Type{T}, parameters::Parameters,enum::IssmEnum) where T #{{{

	param = parameters.lookup[enum]
	return GetParameterValue(param)::T

end#}}}
