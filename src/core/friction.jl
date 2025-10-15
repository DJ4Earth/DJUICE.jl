#Friction class definition

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
struct  CoreWeertmanFriction <: CoreFriction#{{{
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
mutable struct CoreFluxDNNFriction <: CoreFriction#{{{
	dnnChain::Vector{Flux.Chain}
	dtx::Vector{StatsBase.ZScoreTransform}
	dty::Vector{StatsBase.ZScoreTransform}
	xyz_list::Matrix{Float64}
	vx_input::Input
	vy_input::Input
   b_input::Input
   H_input::Input
	rho_ice::Float64
	rho_water::Float64
	g::Float64
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
		dnnChain  = FindParam(Vector{Flux.Chain{}}, element, FrictionDNNChainEnum)
		dtx  = FindParam(Vector{StatsBase.ZScoreTransform{Float64, Vector{Float64}} }, element, FrictionDNNdtxEnum)
		dty  = FindParam(Vector{StatsBase.ZScoreTransform{Float64, Vector{Float64}} }, element, FrictionDNNdtyEnum)
		H_input   = GetInput(element, ThicknessEnum)
		b_input   = GetInput(element, BaseEnum)
		ssx_input   = GetInput(element, SurfaceSlopeXEnum)
		ssy_input   = GetInput(element, SurfaceSlopeYEnum)
		bsx_input   = GetInput(element, BedSlopeXEnum)
		bsy_input   = GetInput(element, BedSlopeYEnum)

      xyz_list = GetVerticesCoordinates(element.vertices)
		
		rho_ice   = FindParam(Float64, element, MaterialsRhoIceEnum)
		rho_water = FindParam(Float64, element, MaterialsRhoSeawaterEnum)
		g         = FindParam(Float64, element, ConstantsGEnum)

		return CoreFluxDNNFriction(dnnChain,dtx,dty,xyz_list,vx_input,vy_input,b_input,H_input,rho_ice,rho_water,g)
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
	
	if vmag==0.0 && m<1.0
		return 0.0
	else
		return c^2*vmag^(m-1)
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
@inline function Alpha2(friction::CoreFluxDNNFriction, gauss::GaussTria, i::Int64)#{{{
	bed = GetInputValue(friction.b_input, gauss, i)
	H = GetInputValue(friction.H_input, gauss, i)
	vx = GetInputValue(friction.vx_input, gauss, i)
	vy = GetInputValue(friction.vy_input, gauss, i)
	h = bed + H

	# Get the velocity
	vmag = VelMag(friction, gauss, i)

	# velocity gradients
	dvx = GetInputDerivativeValue(friction.vx_input,friction.xyz_list,gauss,i)
	dvy = GetInputDerivativeValue(friction.vy_input,friction.xyz_list,gauss,i)
	vxdx = dvx[1]
	vxdy = dvx[2]
	vydx = dvy[1]
	vydy = dvy[2]

	# Get effective pressure
	Neff = EffectivePressure(friction, gauss, i)

	# need to change according to the construction of DNN
	alpha2 = 0.0
	for i in 1:length(friction.dnnChain)
		xin = StatsBase.transform(friction.dtx[i], (reshape(vcat(vx, vy, vxdx, vxdy, vydx, vydy, bed, h), 8, :)))
		pred = StatsBase.reconstruct(friction.dty[i], friction.dnnChain[i](xin))
		alpha2 += first(pred)
	end
	# Average
	alpha2 = alpha2 / length(friction.dnnChain)
	if ( (vmag == 0.0) | (alpha2 < 0.0) )
		alpha2 = 0.0
	else
		alpha2 = alpha2 ./ vmag
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
