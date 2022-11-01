#Friction class definition

abstract type CoreFriction end
struct CoreBuddFriction <: CoreFriction #{{{
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
struct CoreWeertmanFriction <: CoreFriction#{{{
   c_input::Input
	vx_input::Input
	vy_input::Input
	m_input::Input
end# }}}
struct CoreDNNFriction <: CoreFriction#{{{
	dnnChain::Flux.Chain
	dtx::StatsBase.ZScoreTransform
	dty::StatsBase.ZScoreTransform
	vx_input::Input
	vy_input::Input
   c_input::Input
   b_input::Input
   H_input::Input
   ssx_input::Input
   ssy_input::Input
   bsx_input::Input
   bsy_input::Input
	rho_ice::Float64
	rho_water::Float64
	g::Float64
	Cmax::Float64
	velThreshold::Float64
end# }}}

function CoreFriction(element::Tria) #{{{

	frictionlaw = FindParam(Int64, element, FrictionLawEnum)
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
	elseif frictionlaw==10
		dnnChain  = FindParam(Flux.Chain, element, FrictionDNNChainEnum)
		dtx  = FindParam(StatsBase.ZScoreTransform, element, FrictionDNNdtxEnum)
		dty  = FindParam(StatsBase.ZScoreTransform, element, FrictionDNNdtyEnum)
		c_input   = GetInput(element, FrictionCoefficientEnum)
		H_input   = GetInput(element, ThicknessEnum)
		b_input   = GetInput(element, BaseEnum)
		ssx_input   = GetInput(element, SurfaceSlopeXEnum)
		ssy_input   = GetInput(element, SurfaceSlopeYEnum)
		bsx_input   = GetInput(element, BedSlopeXEnum)
		bsy_input   = GetInput(element, BedSlopeYEnum)

		rho_ice   = FindParam(Float64, element, MaterialsRhoIceEnum)
		rho_water = FindParam(Float64, element, MaterialsRhoSeawaterEnum)
		g         = FindParam(Float64, element, ConstantsGEnum)

		Cmax          = FindParam(Float64, element, FrictionCmaxEnum)
		velThreshold = FindParam(Float64, element, VelThresholdEnum)

		return CoreDNNFriction(dnnChain,dtx,dty,vx_input,vy_input,c_input,b_input,H_input,ssx_input,ssy_input,bsx_input,bsy_input, rho_ice, rho_water, g, Cmax, velThreshold)
	else
		error("Friction ",typeof(md.friction)," not supported yet")
	end
end#}}}

#vertices functions
function Alpha2(friction::CoreBuddFriction, gauss::GaussTria, i::Int64) #{{{

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
function Alpha2(friction::CoreWeertmanFriction, gauss::GaussTria, i::Int64)#{{{
	c = GetInputValue(friction.c_input, gauss, i)
	m = GetInputValue(friction.m_input, gauss, i)
	vmag = VelMag(friction, gauss, i)
	
	if vmag==0.0 && m<1.0
		return 0.0
	else
		return c^2*vmag^(m-1)
	end
end#}}}
function Alpha2(friction::CoreDNNFriction, gauss::GaussTria, i::Int64)#{{{
	b = GetInputValue(friction.b_input, gauss, i)
	H = GetInputValue(friction.H_input, gauss, i)
	vx = GetInputValue(friction.vx_input, gauss, i)
	vy = GetInputValue(friction.vy_input, gauss, i)
	ssx = GetInputValue(friction.ssx_input, gauss, i)
	ssy = GetInputValue(friction.ssy_input, gauss, i)
	bsx = GetInputValue(friction.bsx_input, gauss, i)
	bsy = GetInputValue(friction.bsy_input, gauss, i)

	# Get the velocity
	vmag = VelMag(friction, gauss, i)
	# Get effective pressure
	Neff = EffectivePressure(friction, gauss, i)

	# need to change according to the construction of DNN
	xin = StatsBase.transform(friction.dtx, (reshape(vcat(vmag, b, H, ssx, ssy, bsx, bsy), 7, :)))
	pred = StatsBase.reconstruct(friction.dty, friction.dnnChain(xin))
	alpha2 = first(pred)
	if ( (vmag == 0.0) | (alpha2 < 0.0) )
		alpha2 = 0.0
	elseif vmag > friction.velThreshold
		alpha2 = friction.Cmax .* Neff ./ vmag
	else
		alpha2 = alpha2 ./ vmag
	end
	return alpha2
end#}}}
function VelMag(friction::CoreFriction, gauss::GaussTria, i::Int64) #{{{
	vx = GetInputValue(friction.vx_input, gauss, i)
	vy = GetInputValue(friction.vy_input, gauss, i)
	vmag = sqrt(vx^2+vy^2)
end #}}}
function EffectivePressure(friction::CoreFriction, gauss::GaussTria, i::Int64) #{{{
	# Get effective pressure
	H = GetInputValue(friction.H_input, gauss, i)
	b = GetInputValue(friction.b_input, gauss, i)
	N = friction.rho_ice*friction.g*H + friction.rho_water*friction.g*b
	if(N<0.0) N=0.0 end
	return N
end #}}}
