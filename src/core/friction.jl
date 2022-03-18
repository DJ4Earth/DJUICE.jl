#Friction class definition

abstract type CoreFriction end
struct CoreBuddFriction#{{{
	H_input::ElementInput
	b_input::ElementInput
	c_input::ElementInput
	rho_ice::Float64
	rho_water::Float64
	g::Float64
end# }}}
struct CoreWeertmanFriction#{{{
   c_input::ElementInput
	m_input::ElementInput
	vx_input::ElementInput
	vy_input::ElementInput
end# }}}

function CoreFriction(element::Tria) #{{{

	frictionlaw = FindParam(Int64, element, FrictionLawEnum)

	if frictionlaw==1
		H_input  = GetInput(element, ThicknessEnum)
		b_input  = GetInput(element, BaseEnum)
		c_input  = GetInput(element, FrictionCoefficientEnum)

		rho_ice   = FindParam(Float64, element, MaterialsRhoIceEnum)
		rho_water = FindParam(Float64, element, MaterialsRhoSeawaterEnum)
		g         = FindParam(Float64, element, ConstantsGEnum)

		return CoreBuddFriction(H_input, b_input, c_input, rho_ice, rho_water, g)
	elseif frictionlaw==2
		c_input   = GetInput(element, FrictionCEnum)
		m_input   = GetInput(element, FrictionMEnum)
		vx_input  = GetInput(element, VxEnum)
		vy_input  = GetInput(element, VyEnum)
		return CoreWeertmanFriction(c_input,m_input,vx_input,vy_input)
	else
		error("Friction ",typeof(md.friction)," not supported yet")
	end
end#}}}

#vertices functions
function Alpha2(friction::CoreBuddFriction, gauss::GaussTria, i::Int64) #{{{

	#Get effective pressure
	H = GetInputValue(friction.H_input, gauss, i)
	b = GetInputValue(friction.b_input, gauss, i)
	c = GetInputValue(friction.c_input, gauss, i)
	N = friction.rho_ice*friction.g*H + friction.rho_water*friction.g*b

	if(N<0.0) N=0.0 end

	return c^2*N
end #}}}
function Alpha2(friction::CoreWeertmanFriction, gauss::GaussTria, i::Int64)#{{{
	c = GetInputValue(friction.c_input, gauss, i)
	m = GetInputValue(friction.m_input, gauss, i)
	vx = GetInputValue(friction.vx_input, gauss, i)
	vy = GetInputValue(friction.vy_input, gauss, i)
	
	if sqrt(vx^2+vy^2)==0.0 && m<1.0
		return 0.0
	else
		return c^2*sqrt(vx^2+vy^2)^(m-1)
	end
end#}}}
