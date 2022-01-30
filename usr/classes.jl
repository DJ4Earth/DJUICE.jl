using Printf

#Model fields
#Mesh {{{
abstract type AbstractMesh end
mutable struct Mesh2dTriangle <: AbstractMesh
	numberofvertices::Int64
	numberofelements::Int64
	x::Vector{Float64}
	y::Vector{Float64}
	elements::Matrix{Int64}
	segments::Matrix{Int64}
	vertexonboundary::Vector{Bool}
end
function Mesh2dTriangle() #{{{
	return Mesh2dTriangle( 0, 0, Vector{Float64}(undef,0), Vector{Float64}(undef, 0), Matrix{Int64}(undef, 0, 0), Matrix{Int64}(undef, 0, 0), Vector{Bool}(undef,0))
end# }}}
function Base.show(io::IO, this::Mesh2dTriangle)# {{{
	IssmStructDisp(io, this)
end# }}}
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
function Mesh3dPrism() #{{{
	return Mesh3dPrism( 0, 0, 0, Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0), Matrix{Int64}(undef, 0, 0), Matrix{Int64}(undef, 0, 0), Vector{Bool}(undef,0))
end# }}}
#}}}
#Geometry{{{
mutable struct Geometry
	surface::Vector{Float64}
	base::Vector{Float64}
	thickness::Vector{Float64}
	bed::Vector{Float64}
end
function Geometry() #{{{
	return Geometry( Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0), Vector{Float64}(undef,0))
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
#Friction {{{
abstract type AbstractFriction end
mutable struct BuddFriction <: AbstractFriction
	coefficient::Vector{Float64}
end
function BuddFriction() #{{{
	return BuddFriction(Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::BuddFriction)# {{{
	IssmStructDisp(io, this)
end# }}}
mutable struct WeertmanFriction <: AbstractFriction
	C::Vector{Float64}
	m::Vector{Float64}
end
function WeertmanFriction() #{{{
	return WeertmanFriction(Vector{Float64}(undef,0),Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::WeertmanFriction)# {{{
   IssmStructDisp(io, this)
end# }}}
# }}}
#Basalforcings {{{
mutable struct Basalforcings
	groundedice_melting_rate::Vector{Float64}
	floatingice_melting_rate::Vector{Float64}
end
function Basalforcings() #{{{
	return Basalforcings( Vector{Float64}(undef,0), Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::Basalforcings)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Surfaceforcings {{{
mutable struct SMBforcings
	mass_balance::Vector{Float64}
end
function SMBforcings() #{{{
	return SMBforcings( Vector{Float64}(undef,0))
end# }}}
function Base.show(io::IO, this::SMBforcings)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}
#Timestepping{{{
abstract type AbstractTimestepping end
mutable struct Timestepping <: AbstractTimestepping
	start_time::Float64
	final_time::Float64
	time_step::Float64
end
function Timestepping() #{{{
	return Timestepping( 0., 0., 0.)
end# }}}
function Base.show(io::IO, this::Timestepping)# {{{
	IssmStructDisp(io, this)
end# }}}
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
	return Transient( true, true, true, true, true)
end# }}}
function Base.show(io::IO, this::Transient)# {{{
	IssmStructDisp(io, this)
end# }}}
# }}}

#Model structure
mutable struct model
	mesh::AbstractMesh
	geometry::Geometry
	mask::Mask
	materials::Materials
	initialization::Initialization
	stressbalance::Stressbalance
	constants::Constants
	results::Dict
	friction::AbstractFriction
	basalforcings::Basalforcings
	smb::SMBforcings
	timestepping::Timestepping
	masstransport::Masstransport
	transient::Transient
end
function model() #{{{
	return model( Mesh2dTriangle(), Geometry(), Mask(), Materials(),
					 Initialization(),Stressbalance(), Constants(), Dict(),
					 BuddFriction(), Basalforcings(), SMBforcings(), Timestepping(),
					 Masstransport(), Transient())
end#}}}
function model(matmd::Dict) #{{{

	#initialize output
	md = model()

	#Loop over all possible fields
	for name1 in keys(matmd)
		if !(Symbol(name1) in fieldnames(model))
			println("could not recover md.",name1)
			continue
		end
		mdfield  = getfield(md,Symbol(name1))
		matfield = matmd[name1]
		for name2 in keys(matfield)
			if !(Symbol(name2) in fieldnames(typeof(mdfield)))
				println("could not recover md.",name1,".",name2)
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
				error("Don't know how to convert ",typeof(value_matlab)," to ",typeof(value_julia))
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

end# }}}
