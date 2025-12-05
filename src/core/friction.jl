#Friction class definition
using Lux

abstract type CoreFriction end
mutable struct CoreBuddFriction <: CoreFriction #{{{
	c_input::Input
	vx_input::Input
	vy_input::Input
	p_input::Input
	q_input::Input
	H_input::Input
	b_input::Input
	rho_ice::Float64
	rho_water::Float64
	g::Float64
end# }}}
mutable struct CoreWeertmanFriction <: CoreFriction#{{{
   c_input::Input
	vx_input::Input
	vy_input::Input
	m_input::Input
end# }}}
mutable struct CoreSchoofFriction <: CoreFriction #{{{
	c_input::Input
	vx_input::Input
	vy_input::Input
	m_input::Input
	Cmax_input::Input
	H_input::Input
	b_input::Input
	rho_ice::Float64
	rho_water::Float64
	g::Float64
end# }}}
mutable struct CoreDNNFriction <: CoreFriction#{{{
	model::AbstractLuxLayer
	ps
	st
	inputScale::Float64
	outputScale::Float64
   c_input::Input
	vx_input::Input
	vy_input::Input
end# }}}

function CoreFriction(element::Tria, ::Val{frictionlaw}) where frictionlaw #{{{

	vx_input  = GetInput(element, VxEnum)
	vy_input  = GetInput(element, VyEnum)

	if frictionlaw==1
		H_input  = GetInput(element, ThicknessEnum)
		b_input  = GetInput(element, BaseEnum)
		c_input  = GetInput(element, FrictionCoefficientEnum)
		p_input  = GetInput(element, FrictionPEnum)
		q_input  = GetInput(element, FrictionQEnum)

		rho_ice   = FindParam(Float64, element, MaterialsRhoIceEnum)
		rho_water = FindParam(Float64, element, MaterialsRhoSeawaterEnum)
		g         = FindParam(Float64, element, ConstantsGEnum)

		return CoreBuddFriction(c_input, vx_input, vy_input, p_input, q_input, H_input, b_input, rho_ice, rho_water, g)
	elseif frictionlaw==2
		c_input   = GetInput(element, FrictionCEnum)
		m_input   = GetInput(element, FrictionMEnum)
		return CoreWeertmanFriction(c_input,vx_input,vy_input,m_input)
	elseif frictionlaw==11
		H_input  = GetInput(element, ThicknessEnum)
		b_input  = GetInput(element, BaseEnum)
		c_input   = GetInput(element, FrictionCEnum)
		m_input   = GetInput(element, FrictionMEnum)
		Cmax_input   = GetInput(element, FrictionCmaxEnum)

		rho_ice   = FindParam(Float64, element, MaterialsRhoIceEnum)
		rho_water = FindParam(Float64, element, MaterialsRhoSeawaterEnum)
		g         = FindParam(Float64, element, ConstantsGEnum)

		return CoreSchoofFriction(c_input, vx_input, vy_input, m_input, Cmax_input, H_input, b_input, rho_ice, rho_water, g)
	elseif frictionlaw==20
		c_input   = GetInput(element, FrictionCEnum)
		model  = FindParam(AbstractLuxLayer, element, FrictionDNNEnum)
		ps  = FindParam(NamedTuple, element, FrictionDNNpsEnum)
		st  = FindParam(NamedTuple, element, FrictionDNNstEnum)
		inputScale  = FindParam(Float64, element, FrictionDNNInputScaleEnum)
		outputScale  = FindParam(Float64, element, FrictionDNNOutputScaleEnum)

		return CoreDNNFriction(model,ps,st,inputScale,outputScale,c_input,vx_input,vy_input)
	else
		error("Friction ",typeof(md.friction)," not supported yet")
	end
end#}}}

#vertices functions
@inline function Alpha2(friction::CoreBuddFriction, gauss::GaussTria, i::Int64) #{{{

	# Recover parameters
	p = GetInputValue(friction.p_input, gauss, i)
	q = GetInputValue(friction.q_input, gauss, i)

	# Compute r and s coefficients
	r = q / p
	s = 1.0/p
	c = GetInputValue(friction.c_input, gauss, i)
	
	# Get effective pressure
	N = EffectivePressure(friction, gauss, i)
	# Get the velocity
	vmag = VelMag(friction, gauss, i)

	if(N<0.0) N=0.0 end

	if (s == 1.0)
		alpha2 = c^2 * (N^r)
	else
		if (vmag == 0.0 && s<1.0)
			alpha2 = 0.0
		else
			alpha2 = c^2 * (N^r) * vmag^(s-1.0)
		end
	end
	return alpha2
end #}}}
@inline function Alpha2(friction::CoreWeertmanFriction, gauss::GaussTria, i::Int64)#{{{
	c = GetInputValue(friction.c_input, gauss, i)
	m = GetInputValue(friction.m_input, gauss, i)
	vmag = VelMag(friction, gauss, i)
	
	if vmag==0.0 && (1.0/m)<1.0
		return 0.0
	else
		return c^2*vmag^(1.0/m-1.0)
	end
end#}}}
@inline function Alpha2(friction::CoreSchoofFriction, gauss::GaussTria, i::Int64) #{{{

	# Recover parameters
	m = GetInputValue(friction.m_input, gauss, i)
	c = GetInputValue(friction.c_input, gauss, i)
	Cmax = GetInputValue(friction.Cmax_input, gauss, i)

	# Get effective pressure
	N = EffectivePressure(friction, gauss, i)
	# Get the velocity
	vmag = VelMag(friction, gauss, i)

	if ((vmag<1.0e-20) || (N == 0.0))
		alpha2 = 0.0
	else
		alpha2 = (c^2 * vmag^(m-1)) / ((1.0 + (c^2/(Cmax*N))^(1.0/m)*vmag)^m)
	end
	return alpha2
end #}}}
@inline function Alpha2(friction::CoreDNNFriction, gauss::GaussTria, i::Int64)#{{{
	vx = GetInputValue(friction.vx_input, gauss, i)
	vy = GetInputValue(friction.vy_input, gauss, i)

	c = GetInputValue(friction.c_input, gauss, i)
	# Get the velocity
	vmag = VelMag(friction, gauss, i)
	# construct nn model
	smodel = StatefulLuxLayer(friction.model, friction.ps, friction.st)

	if (vmag == 0.0 )
		alpha2 = 0.0
	else
		alpha2 = c^2*((smodel([vmag]/friction.inputScale)[1])*friction.outputScale)/vmag
	end
	return alpha2
end#}}}
@inline function VelMag(friction::CoreFriction, gauss::GaussTria, i::Int64) #{{{
	vx = GetInputValue(friction.vx_input, gauss, i)
	vy = GetInputValue(friction.vy_input, gauss, i)
	vmag = sqrt(vx^2+vy^2)
end #}}}
@inline function EffectivePressure(friction::CoreFriction, gauss::GaussTria, i::Int64) #{{{
	# Get effective pressure
	H = GetInputValue(friction.H_input, gauss, i)
	b = GetInputValue(friction.b_input, gauss, i)
	N = friction.rho_ice*friction.g*H + friction.rho_water*friction.g*b
	if(N<0.0) N=0.0 end
	return N
end #}}}
