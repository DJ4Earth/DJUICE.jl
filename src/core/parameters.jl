#Parameter class definition
abstract type Parameter end
struct BoolParam <: Parameter #{{{
	enum::IssmEnum
	value::Bool
end# }}}
struct DoubleParam <: Parameter #{{{
	enum::IssmEnum
	value::Float64
end# }}}
struct IntParam <: Parameter #{{{
	enum::IssmEnum
	value::Int64
end# }}}
struct StringParam <: Parameter #{{{
	enum::IssmEnum
	value::String
end# }}}
struct SymbolParam <: Parameter #{{{
	enum::IssmEnum
	value::Symbol
end# }}}
mutable struct StringArrayParam <: Parameter #{{{
	enum::IssmEnum
	value::Vector{String}
end# }}}
mutable struct FluxChainParam <: Parameter #{{{
	enum::IssmEnum
	value::Vector{Flux.Chain{}}
end# }}}
mutable struct StatsBaseTransformParam <: Parameter #{{{
	enum::IssmEnum
	value::Vector{StatsBase.ZScoreTransform{Float64, Vector{Float64}} } 
end# }}}

#Parameters dataset class definition
mutable struct Parameters #{{{
	lookup::Dict{IssmEnum,Parameter}
	double_lookup::Dict{IssmEnum,DoubleParam}
end# }}}
function Parameters() #{{{
	return Parameters(Dict{IssmEnum,Parameter}(), Dict{IssmEnum,DoubleParam}())
end # }}}

#Parameter functions
function GetParameterValue(param::BoolParam) #{{{
	return param.value::Bool
end#}}}
function GetParameterValue(param::DoubleParam) #{{{
	return param.value::Float64
end#}}}
function GetParameterValue(param::IntParam) #{{{
	return param.value::Int64
end#}}}
function GetParameterValue(param::StringParam) #{{{
	return param.value::String
end#}}}
function GetParameterValue(param::SymbolParam) #{{{
	return param.value::Symbol
end#}}}
function GetParameterValue(param::StringArrayParam) #{{{
	return param.value::Vector{String}
end#}}}
function GetParameterValue(param::FluxChainParam) #{{{
	return param.value::Vector{Flux.Chain{}}
end#}}}
function GetParameterValue(param::StatsBaseTransformParam) #{{{
	return param.value::Vector{StatsBase.ZScoreTransform{Float64, Vector{Float64} }}
end#}}}

#Parameters functions
function AddParam(parameters::Parameters,value::Bool, enum::IssmEnum) #{{{

	parameters.lookup[enum] = BoolParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Float64,enum::IssmEnum) #{{{

	parameters.double_lookup[enum] = DoubleParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Int64, enum::IssmEnum) #{{{

	parameters.lookup[enum] = IntParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::String, enum::IssmEnum) #{{{

	parameters.lookup[enum] = StringParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Symbol, enum::IssmEnum) #{{{

	parameters.lookup[enum] = SymbolParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Vector{String}, enum::IssmEnum) #{{{

	parameters.lookup[enum] = StringArrayParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Vector{Flux.Chain{}}, enum::IssmEnum) #{{{

	parameters.lookup[enum] = FluxChainParam(enum,value)

	return nothing
end#}}}
function AddParam(parameters::Parameters,value::Vector{StatsBase.ZScoreTransform{Float64, Vector{Float64}} }, enum::IssmEnum) #{{{

	parameters.lookup[enum] = StatsBaseTransformParam(enum,value)

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
