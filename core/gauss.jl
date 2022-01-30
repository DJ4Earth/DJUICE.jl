#Gauss class definition
struct GaussTria #{{{
	numgauss::Int64
	weights::Vector{Float64}
	coords1::Vector{Float64}
	coords2::Vector{Float64}
	coords3::Vector{Float64}
end #}}}
function Base.show(io::IO, this::GaussTria)# {{{

	println(io,"GaussTria:")
	println(io,"   numgauss: ",this.numgauss)
	println(io,"   weights:  ",this.weights)
	println(io,"   coords1:  ",this.coords1)
	println(io,"   coords2:  ",this.coords2)
	println(io,"   coords3:  ",this.coords3)
end# }}}

#Gauss constructor
function GaussTria(order::Int64) #{{{

	#=Gauss quadrature points for the triangle.
	Higher-order points from D.A. Dunavant, "High Degree Efficient
	Symmetrical Gaussian Quadrature Rules for the Triangle", IJNME,
	Vol. 21, pp. 1129-1148 (1985), as transcribed for Probe rules3.=#

	if(order==1)
		npoints = 1
		weights = [1.732050807568877]
		coords1 = [0.333333333333333]
		coords2 = [0.333333333333333]
		coords3 = [0.333333333333333]
	elseif(order==2)
		npoints = 3
		weights = [0.577350269189625; 0.577350269189625; 0.577350269189625]
		coords1 = [0.666666666666667; 0.166666666666667; 0.166666666666667]
		coords2 = [0.166666666666667; 0.666666666666667; 0.166666666666667]
		coords3 = [0.166666666666667; 0.166666666666667; 0.666666666666667]
	elseif(order==3)
		npoints = 4
		weights = [-0.974278579257493; 0.902109795608790; 0.902109795608790; 0.902109795608790]
		coords1 = [ 0.333333333333333; 0.600000000000000; 0.200000000000000; 0.200000000000000]
		coords2 = [ 0.333333333333333; 0.200000000000000; 0.600000000000000; 0.200000000000000]
		coords3 = [ 0.333333333333333; 0.200000000000000; 0.200000000000000; 0.600000000000000]
	else
		error("order ",order," not supported yet");
	end

	return GaussTria(npoints,weights,coords1,coords2,coords3)
end# }}}
function GaussTria(finiteelement::IssmEnum) #{{{

	if(finiteelement==P0Enum)
		npoints = 1
		weights = [1.]
		coords1 = [0.333333333333333]
		coords2 = [0.333333333333333]
		coords3 = [0.333333333333333]
	elseif(finiteelement==P1Enum)
			npoints = 3
			weights = 0.333333333333333*ones(3)
			coords1 = [1.; 0.; 0.]
			coords2 = [0.; 1.; 0.]
			coords3 = [0.; 0.; 1.]
	else
		error("finite element ", finiteelement," not supported yet");
	end

	return GaussTria(npoints,weights,coords1,coords2,coords3)
end# }}}
function GaussTria(area_coordinates::Matrix{Float64}, order::Int64) #{{{
	#=Gauss-Legendre quadrature points.

	The recurrence coefficients for Legendre polynomials on (-1,1)
	are defined (from the ORTHPOL subroutine RECUR with ipoly=1) as:

	alpha(i)=0.
	beta (i)=1./(4.-1./(i-1)^2))

	For degree p, the required number of Gauss-Legendre points is
	n>=(p+1)/2.=#

	if(order==1)
		npoint  = 1
		weights = [2.000000000000000]
		coords  = [0.000000000000000]
	elseif(order==2)
		npoints = 2
		weights = [1.000000000000000, 1.000000000000000]
		coords  = [-0.577350269189626, 0.577350269189626]
	elseif(order==3)
		npoints = 3
		weights = [0.555555555555556, 0.888888888888889, 0.555555555555556]
		coords  = [-0.774596669241483, 0.000000000000000, 0.774596669241483]
	elseif(order==4)
		npoints = 4
		weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]
		coords  = [-0.861136311594053,-0.339981043584856, 0.339981043584856, 0.861136311594053]
	else
      error("order ",order," not supported yet");
	end

   coords1  = Vector{Float64}(undef,npoints)
   coords2  = Vector{Float64}(undef,npoints)
   coords3  = Vector{Float64}(undef,npoints)
   for i in 1:npoints
      coords1[i]  = 0.5*(area_coordinates[1,1]+area_coordinates[2,1]) + 0.5*coords[i]*(area_coordinates[2,1]-area_coordinates[1,1]);
      coords2[i]  = 0.5*(area_coordinates[1,2]+area_coordinates[2,2]) + 0.5*coords[i]*(area_coordinates[2,2]-area_coordinates[1,2]);
      coords3[i]  = 0.5*(area_coordinates[1,3]+area_coordinates[2,3]) + 0.5*coords[i]*(area_coordinates[2,3]-area_coordinates[1,3]);
   end

	return GaussTria(npoints, weights, coords1, coords2, coords3)
end# }}}
