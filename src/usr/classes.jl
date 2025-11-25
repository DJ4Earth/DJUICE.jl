using Printf
using Lux
using StatsBase

#Model fields
#Mesh (Abstract){{{
abstract type AbstractMesh end
function Base.show(io::IO, this::AbstractMesh) IssmStructDisp(io, this) end

#Mesh2dTriangle
mutable struct Mesh2dTriangle <: AbstractMesh
	numberofvertices::Int64
	numberofelements::Int64
	x::Vector{Float64}
	y::Vector{Float64}
	elements::Matrix{Int64}
	segments::Matrix{Int64}
	vertexonboundary::Vector{Bool}
end
function Mesh2dTriangle()
	return Mesh2dTriangle( 0, 0, Vector{Float64}(undef,0), Vector{Float64}(undef, 0), Matrix{Int64}(undef, 0, 0), Matrix{Int64}(undef, 0, 0), Vector{Bool}(undef,0))
end

#Mesh3dPrism
mutable struct Mesh3dPrism{T} <: AbstractMesh
	numberofvertices::Int64
	numberofelements::Int64
	numberoflayers::Int64
	x::Vector{Float64}
	y::Vector{Float64}
	z::Vector{Float64}
	elements::Matrix{Int64}
	segments::Matrix{Int64}
	vertexonboundary::Vector{Bool}
end
function Mesh3dPrism()
	return Mesh3dPrism(0, 0, 0,
							 Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0),
							 Matrix{Int64}(undef, 0, 0), Matrix{Int64}(undef, 0, 0), Vector{Bool}(undef,0))
end
#}}}
#Geometry{{{
mutable struct Geometry
	surface::Vector{Float64}
	base::Vector{Float64}
	thickness::Vector{Float64}
	bed::Vector{Float64}
	ssx::Vector{Float64}
	ssy::Vector{Float64}
	bsx::Vector{Float64}
	bsy::Vector{Float64}
end
function Geometry() #{{{
	return Geometry( Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0), 
						 Vector{Float64}(undef,0), Vector{Float64}(undef,0),
						 Vector{Float64}(undef,0), Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::Geometry)# {{{
	IssmStructDisp(io, this)
end# }}}
#}}}
#Mask {{{
mutable struct Mask
	ocean_levelset::Vector{Float64}
	ice_levelset::Vector{Float64}
end
function Mask() #{{{
	return Mask( Vector{Float64}(undef,0), Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::Mask)# {{{
	IssmStructDisp(io, this)
end# }}}
#}}}
#Initialization{{{
mutable struct Initialization
	vx::Vector{Float64}
	vy::Vector{Float64}
end
function Initialization() #{{{
	return Initialization( Vector{Float64}(undef,0), Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::Initialization)# {{{
	IssmStructDisp(io, this)
end# }}}
#}}}
#Stressbalance {{{
mutable struct Stressbalance
	spcvx::Vector{Float64}
	spcvy::Vector{Float64}
	restol::Float64
	reltol::Float64
	abstol::Float64
	maxiter::Int64
end
function Stressbalance() #{{{
	return Stressbalance( Vector{Float64}(undef,0), Vector{Float64}(undef,0), 1.e-4, 0.01, 10., 100)
end# }}}
function Base.show(io::IO, this::Stressbalance)# {{{
	IssmStructDisp(io, this)
end# }}}
#}}}
#Constants{{{
mutable struct Constants
	g::Float64
	yts::Float64
end
function Constants() #{{{
	return Constants( 9.81,  365*24*3600.)
end# }}}
function Base.show(io::IO, this::Constants)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Materials {{{
mutable struct Materials
	rho_ice::Float64
	rho_water::Float64
	rho_freshwater::Float64
	mu_water::Float64
	heatcapacity::Float64
	latentheat::Float64
	thermalconductivity::Float64
	temperateiceconductivity::Float64
	effectiveconductivity_averaging::Int64
	meltingpoint::Float64
	beta::Float64
	mixed_layer_capacity::Float64
	thermal_exchange_velocity::Float64
	rheology_B::Vector{Float64}
	rheology_n::Vector{Float64}
	rheology_law::String
end
function Materials() #{{{
	return Materials(917., 1023., 1000., 0.001787, 2093., 3.34*10^5, 2.4, .24, 1, 273.15, 9.8*10^-8, 3974., 1.00*10^-4, Vector{Float64}(undef,0), Vector{Float64}(undef,0), "Cuffey")
end# }}}
function Base.show(io::IO, this::Materials)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Friction (Abstract){{{
abstract type AbstractFriction end
function Base.show(io::IO, this::AbstractFriction) IssmStructDisp(io, this) end

#BuddFriction
mutable struct BuddFriction <: AbstractFriction
	coefficient::Vector{Float64}
	p::Vector{Float64}
	q::Vector{Float64}
end
function BuddFriction()
	return BuddFriction(Vector{Float64}(undef,0),Vector{Float64}(undef,0),Vector{Float64}(undef,0))
end

#WeertmanFriction
mutable struct WeertmanFriction <: AbstractFriction
	C::Vector{Float64}
	m::Vector{Float64}
end
function WeertmanFriction()
	return WeertmanFriction(Vector{Float64}(undef,0),Vector{Float64}(undef,0))
end

#SchoofFriction
mutable struct SchoofFriction <: AbstractFriction
	C::Vector{Float64}
	m::Vector{Float64}
	Cmax::Vector{Float64}
end
function SchoofFriction()
	return SchoofFriction(Vector{Float64}(undef,0),Vector{Float64}(undef,0),Vector{Float64}(undef,0))
end

#DNNFriction: using lux
mutable struct DNNFriction <: AbstractFriction
	C::Vector{Float64}
	model::AbstractLuxLayer
	ps
	st
	input_scale::Float64
	output_scale::Float64
end
function DNNFriction() 
	return DNNFriction(Vector{Float64}(undef,0),
							 Lux.Chain(),
							 NamedTuple{},
							 NamedTuple{},
							 1.0, 1.0)
						 end
						 # }}}
#Basalforcings (Abstract) {{{
abstract type AbstractBasalforcings end
function Base.show(io::IO, this::AbstractBasalforcings) IssmStructDisp(io, this) end

#DefaultBasalforcings
mutable struct DefaultBasalforcings  <: AbstractBasalforcings
	groundedice_melting_rate::Vector{Float64}
	floatingice_melting_rate::Vector{Float64}
end
function DefaultBasalforcings()
	return DefaultBasalforcings( Vector{Float64}(undef,0), Vector{Float64}(undef,0))
end

#LinearBasalforcings
mutable struct LinearBasalforcings  <: AbstractBasalforcings
	deepwater_melting_rate::Float64
	upperwater_melting_rate::Float64
	deepwater_elevation::Float64
	upperwater_elevation::Float64
	groundedice_melting_rate::Vector{Float64}
	perturbation_melting_rate::Vector{Float64}
	geothermalflux::Vector{Float64}
end
function LinearBasalforcings()
	return LinearBasalforcings(50., 0., -800., -400., Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0))
end
# }}}
#Surfaceforcings {{{
mutable struct SMBforcings
	mass_balance::Union{Vector{Float64},Matrix{Float64}}
end
function SMBforcings() #{{{
	return SMBforcings(Vector{Float64}(undef, 0))
end# }}}
function Base.show(io::IO, this::SMBforcings)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Timestepping (Abstract){{{
abstract type AbstractTimestepping end
function Base.show(io::IO, this::AbstractTimestepping) IssmStructDisp(io, this) end

#DefaultTimestepping
mutable struct DefaultTimestepping <: AbstractTimestepping
	start_time::Float64
	final_time::Float64
	time_step::Float64
end
function DefaultTimestepping() 
	return DefaultTimestepping( 0., 0., 0.)
end
# }}}
#Masstransport {{{
mutable struct Masstransport
	spcthickness::Vector{Float64}
	min_thickness::Float64
	stabilization::Int64
end
function Masstransport() #{{{
	return Masstransport( Vector{Float64}(undef,0), 10.0, 1)
end# }}}
function Base.show(io::IO, this::Masstransport)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Transient {{{
mutable struct Transient
	issmb::Bool
	ismasstransport::Bool
	isstressbalance::Bool
	isgroundingline::Bool
	ismovingfront::Bool
end
function Transient() #{{{
	return Transient( true, true, true, false, false)
end# }}}
function Base.show(io::IO, this::Transient)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Inversion{{{
mutable struct Inversion
	iscontrol::Bool
	onlygrad::Bool
	vx_obs::Vector{Float64}
	vy_obs::Vector{Float64}
	min_parameters::Vector{Float64}
	max_parameters::Vector{Float64}
	independent::Vector{Float64}
	maxiter::Int64
	tol::Float64
	independent_string::String
	dependent_string::Vector{String}
end
function Inversion() #{{{
	return Inversion( false, true, Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0), 0, 0., "Friction", Vector{String}(undef,0))
end# }}}
function Base.show(io::IO, this::Inversion)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Calving (Abstract){{{
abstract type AbstractCalving end
function Base.show(io::IO, this::AbstractCalving) IssmStructDisp(io, this) end

#DefaultCalving
mutable struct DefaultCalving <: AbstractCalving
	calvingrate::Vector{Float64}
end
function DefaultCalving()
	return DefaultCalving(Vector{Float64}(undef,0))
end
# }}}
#Levelset{{{
mutable struct Levelset
	spclevelset::Union{Vector{Float64},Matrix{Float64}}
	stabilization::Int64
	reinit_frequency::Int64
	kill_icebergs::Int64
	migration_max::Float64
end
function Levelset() #{{{
	return Levelset(Vector{Float64}(undef,0), 1, 10, 1, 1.0e12)
end# }}}
function Base.show(io::IO, this::Levelset)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Frontalforcings{{{
mutable struct Frontalforcings
	meltingrate::Vector{Float64}
	ablationrate::Vector{Float64}
end
function Frontalforcings() #{{{
	return Frontalforcings(Vector{Float64}(undef,0), Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::Frontalforcings)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Groundingline{{{
mutable struct Groundingline
	migration::String
end
function Groundingline() #{{{
	return Groundingline("None")
end# }}}
function Base.show(io::IO, this::Groundingline)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Verbose{{{
mutable struct Verbose
	mprocessor::Bool
	modules::Bool
	solution::Bool
	solver::Bool
	convergence::Bool
	control::Bool
	autodiff::Bool
end
function Verbose() #{{{
	return Verbose( false, false, false, false, false, false, false)
end# }}}
function Base.show(io::IO, this::Verbose)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}

#Model structure
mutable struct model{Mesh<:AbstractMesh, Friction<:AbstractFriction, Basalforcings<:AbstractBasalforcings, Calving<:AbstractCalving}
	mesh::Mesh
	geometry::Geometry
	mask::Mask
	materials::Materials
	initialization::Initialization
	stressbalance::Stressbalance
	constants::Constants
	results::Dict
	friction::Friction
	basalforcings::Basalforcings
	smb::SMBforcings
	timestepping::DefaultTimestepping
	masstransport::Masstransport
	transient::Transient
	inversion::Inversion
	calving::Calving
	levelset::Levelset
	frontalforcings::Frontalforcings
	groundingline::Groundingline
	verbose::Verbose
end
function model() #{{{
      return model( Mesh2dTriangle(), Geometry(), Mask(), Materials(),
                                       Initialization(),Stressbalance(), Constants(), Dict(),
                                       BuddFriction(), DefaultBasalforcings(), SMBforcings(), DefaultTimestepping(),
													Masstransport(), Transient(), Inversion(), DefaultCalving(), 
													Levelset(), Frontalforcings(), Groundingline(), Verbose())
end#}}}
function model(md::model; mesh::AbstractMesh=md.mesh, friction::AbstractFriction=md.friction, calving::AbstractCalving=md.calving, basalforcings::AbstractBasalforcings=md.basalforcings) #{{{
	return model(mesh, md.geometry, md.mask, md.materials, 
					 md.initialization, md.stressbalance, md.constants, md.results, 
					 friction, basalforcings, md.smb, md.timestepping, 
					 md.masstransport, md.transient, md.inversion, md.calving, 
					 md.levelset, md.frontalforcings, md.groundingline, md.verbose)
end#}}}
function model(matmd::Dict; verbose::Bool=false, friction::AbstractFriction=BuddFriction(), basalforcings::AbstractBasalforcings=DefaultBasalforcings()) #{{{

	#initialize output
	md = model(model(), friction=friction, basalforcings=basalforcings)

	#Loop over all possible fields
	for name1 in keys(matmd)
		if !(Symbol(name1) in fieldnames(model))
			if verbose; println("could not recover md.",name1) end
			continue
		end
		mdfield  = getfield(md,Symbol(name1))
		matfield = matmd[name1]
		for name2 in keys(matfield)
			if !(Symbol(name2) in fieldnames(typeof(mdfield)))
				if verbose; println("could not recover md.",name1,".",name2) end
				continue
			end
			value_matlab = matfield[name2]
			value_julia  = getfield(mdfield, Symbol(name2))

			if typeof(value_matlab)==typeof(value_julia)
				setfield!(mdfield, Symbol(name2), value_matlab)

			elseif typeof(value_matlab)==Float64 && typeof(value_julia)==Int64
				setfield!(mdfield, Symbol(name2), Int64(value_matlab))

			elseif typeof(value_matlab)==Float64 && typeof(value_julia)==Bool
				setfield!(mdfield, Symbol(name2), Bool(value_matlab))

				# TODO: temporarily fix the issue when loading one value from matlab to a vector in Julia
			elseif typeof(value_matlab)==Float64 && typeof(value_julia)==Vector{Float64}
				setfield!(mdfield, Symbol(name2), [value_matlab])

			elseif typeof(value_matlab)==Matrix{Float64} && typeof(value_julia)==Vector{Float64}
				if(size(value_matlab,2)!=1) error("only one column expected") end
				setfield!(mdfield, Symbol(name2), value_matlab[:,1])

			elseif typeof(value_matlab)==Matrix{Float64} && typeof(value_julia)==Matrix{Int64}
				matrix = Matrix{Int64}(undef,size(value_matlab))
				for i in 1:length(value_matlab) matrix[i] = Int64(value_matlab[i]) end
				setfield!(mdfield, Symbol(name2), matrix)

			elseif typeof(value_matlab)==Matrix{Float64} && typeof(value_julia)==Vector{Bool}
				if(size(value_matlab,2)!=1) error("only one column expected") end
				vector = Vector{Bool}(undef,size(value_matlab,1))
				for i in 1:length(vector) vector[i] = Bool(value_matlab[i]) end
				setfield!(mdfield, Symbol(name2), vector)

			else
				error("Don't know how to convert ", name2, " from ",typeof(value_matlab)," to ",typeof(value_julia))
			end
		end
	end

	return md
end#}}}
function Base.show(io::IO, md::model)# {{{

	compact = get(io, :compact, false)

	println(io,"Model:")
	@printf "%19s: %-26s -- %s\n" "mesh" typeof(md.mesh) "mesh properties"
	@printf "%19s: %-26s -- %s\n" "geometry" typeof(md.geometry) "surface elevation, bedrock topography, ice thickness,..."
	@printf "%19s: %-26s -- %s\n" "mask" typeof(md.mask) "defines grounded and floating regions"
	@printf "%19s: %-26s -- %s\n" "materials" typeof(md.materials) "material properties"
	@printf "%19s: %-26s -- %s\n" "initialization" typeof(md.initialization) "initial state"
	@printf "%19s: %-26s -- %s\n" "constants" typeof(md.constants) "physical constants"
	@printf "%19s: %-26s -- %s\n" "friction" typeof(md.friction) "basal friction"
	@printf "%19s: %-26s -- %s\n" "basalforcings" typeof(md.basalforcings) "basal forcings"
	@printf "%19s: %-26s -- %s\n" "smb" typeof(md.smb) "surface mass balance"
	@printf "%19s: %-26s -- %s\n" "timestepping" typeof(md.timestepping) "time stepping for transient simulations"
	@printf "%19s: %-26s -- %s\n" "stressbalance" typeof(md.stressbalance) "parameters stress balance simulations"
	@printf "%19s: %-26s -- %s\n" "masstransport" typeof(md.masstransport) "parameters mass transport simulations"
	@printf "%19s: %-26s -- %s\n" "transient" typeof(md.transient) "parameters for transient simulations"
	@printf "%19s: %-26s -- %s\n" "inversion" typeof(md.inversion) "parameters for inverse methods"
	@printf "%19s: %-26s -- %s\n" "calving" typeof(md.calving) "parameters for calving"
	@printf "%19s: %-26s -- %s\n" "levelset" typeof(md.levelset) "parameters for moving boundaries (level-set method)"
	@printf "%19s: %-26s -- %s\n" "frontalforcings" typeof(md.frontalforcings) "parameters for frontalforcings"
	@printf "%19s: %-26s -- %s\n" "groundingline" typeof(md.groundingline) "parameters for groundingline"
	@printf "%19s: %-26s -- %s\n" "verbose" typeof(md.verbose) "verbosity level in solve"
	@printf "%19s: %-26s -- %s\n" "results" typeof(md.results) "model results"

end# }}}
