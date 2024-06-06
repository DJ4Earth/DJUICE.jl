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

	return nothing
end# }}}

#Gauss constructor
function GaussTria(order::Int64) #{{{
	npoints = GaussLegendreTriaNpoints(Val(order))
	weights = Vector{Float64}(undef, npoints)
	coords1 = Vector{Float64}(undef, npoints)
	coords2 = Vector{Float64}(undef, npoints)
	coords3 = Vector{Float64}(undef, npoints)
	GaussLegendreTria(weights, coords1, coords2, coords3, Val(order))
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
function GaussTria(index::Int64, r1::Float64, r2::Float64, mainlyfloating::Bool, order::Int64) #{{{
	xy_list = Matrix{Float64}(undef,3,2)

	if (mainlyfloating) 
		# Get gauss points
		numgauss = GaussLegendreTriaNpoints(Val(order))
		weights = Vector{Float64}(undef, numgauss)
		coords1 = Vector{Float64}(undef, numgauss)
		coords2 = Vector{Float64}(undef, numgauss)
		coords3 = Vector{Float64}(undef, numgauss)
		GaussLegendreTria(weights, coords1, coords2, coords3, Val(order))

		xy_list[1,1]=0.;  xy_list[1,2]=0.;
		xy_list[2,1]=r1; xy_list[2,2]=0.;
		xy_list[3,1]=0.;  xy_list[3,2]=r2;

		for ii in 1:numgauss
			x = coords1[ii]*xy_list[1,1] + coords2[ii]*xy_list[2,1] + coords3[ii]*xy_list[3,1];
			y = coords1[ii]*xy_list[1,2] + coords2[ii]*xy_list[2,2] + coords3[ii]*xy_list[3,2];

			if (index==1)
				coords1[ii] = 1.0-x-y
				coords2[ii] = x
				coords3[ii] = y
			elseif (index==2)
				coords1[ii] = y
				coords2[ii] = 1.0-x-y
				coords3[ii] = x
			elseif (index==3)
				coords1[ii] = x
				coords2[ii] = y
				coords3[ii] = 1.0-x-y
			else
				error("index ", index, " not supported yet")
			end
			weights[ii] = weights[ii]*r1*r2
		end
	else
		# Double number of gauss points
		gauss1 = GaussTria(order)
		xy_list[1,1]=r1;  xy_list[1,2]=0.;
		xy_list[2,1]=0.; xy_list[2,2]=1.;
		xy_list[3,1]=0.;  xy_list[3,2]=r2;

		for ii in 1:gauss1.numgauss
			x = gauss1.coords1[ii]*xy_list[1,1] + gauss1.coords2[ii]*xy_list[2,1] + gauss1.coords3[ii]*xy_list[3,1];
			y = gauss1.coords1[ii]*xy_list[1,2] + gauss1.coords2[ii]*xy_list[2,2] + gauss1.coords3[ii]*xy_list[3,2];

			if (index==1)
				gauss1.coords1[ii] = 1.0-x-y
				gauss1.coords2[ii] = x
				gauss1.coords3[ii] = y
			elseif (index==2)
				gauss1.coords1[ii] = y
				gauss1.coords2[ii] = 1.0-x-y
				gauss1.coords3[ii] = x
			elseif (index==3)
				gauss1.coords1[ii] = x
				gauss1.coords2[ii] = y
				gauss1.coords3[ii] = 1.0-x-y
			else
				error("index ", index, " not supported yet")
			end
			gauss1.weights[ii] = gauss1.weights[ii]*r1*(1.0-r2)
		end

		gauss2 = GaussTria(order)
		xy_list[1,1]=r1;  xy_list[1,2]=0.;
		xy_list[2,1]=1.; xy_list[2,2]=0.;
		xy_list[3,1]=0.;  xy_list[3,2]=1.;

		for ii in 1:gauss2.numgauss
			x = gauss2.coords1[ii]*xy_list[1,1] + gauss2.coords2[ii]*xy_list[2,1] + gauss2.coords3[ii]*xy_list[3,1];
			y = gauss2.coords1[ii]*xy_list[1,2] + gauss2.coords2[ii]*xy_list[2,2] + gauss2.coords3[ii]*xy_list[3,2];

			if (index==1)
				gauss2.coords1[ii] = 1.0-x-y
				gauss2.coords2[ii] = x
				gauss2.coords3[ii] = y
			elseif (index==2)
				gauss2.coords1[ii] = y
				gauss2.coords2[ii] = 1.0-x-y
				gauss2.coords3[ii] = x
			elseif (index==3)
				gauss2.coords1[ii] = x
				gauss2.coords2[ii] = y
				gauss2.coords3[ii] = 1.0-x-y
			else
				error("index ", index, " not supported yet")
			end
			gauss2.weights[ii] = gauss2.weights[ii]*(1.0-r1)
		end

		numgauss = gauss1.numgauss + gauss2.numgauss
		weights = vcat(gauss1.weights, gauss2.weights)
		coords1 = vcat(gauss1.coords1, gauss2.coords1)
		coords2 = vcat(gauss1.coords2, gauss2.coords2)
		coords3 = vcat(gauss1.coords3, gauss2.coords3)
	end

	return GaussTria(numgauss, weights, coords1, coords2, coords3)
end# }}}

#Numerics
function GaussLegendreTriaNpoints(::Val{order})  where order#{{{
	#=Gauss quadrature points for the triangle.
	Higher-order points from D.A. Dunavant, "High Degree Efficient
	Symmetrical Gaussian Quadrature Rules for the Triangle", IJNME,
	Vol. 21, pp. 1129-1148 (1985), as transcribed for Probe rules3.=#

	if(order==1)
		npoints = 1
	elseif(order==2)
		npoints = 3
	elseif(order==3)
		npoints = 4
	elseif(order==4)
		npoints = 6
	else
		error("order ",order," not supported yet");
	end

	return npoints
end# }}}
function GaussLegendreTria(weights::Vector{Float64}, coords1::Vector{Float64}, coords2::Vector{Float64}, coords3::Vector{Float64}, ::Val{order})  where order#{{{

	#=Gauss quadrature points for the triangle.
	Higher-order points from D.A. Dunavant, "High Degree Efficient
	Symmetrical Gaussian Quadrature Rules for the Triangle", IJNME,
	Vol. 21, pp. 1129-1148 (1985), as transcribed for Probe rules3.=#

	if(order==1)
		npoints = 1
		weights .= [1.732050807568877]
		coords1 .= [0.333333333333333]
		coords2 .= [0.333333333333333]
		coords3 .= [0.333333333333333]
	elseif(order==2)
		npoints = 3
		weights .= [0.577350269189625; 0.577350269189625; 0.577350269189625]
		coords1 .= [0.666666666666667; 0.166666666666667; 0.166666666666667]
		coords2 .= [0.166666666666667; 0.666666666666667; 0.166666666666667]
		coords3 .= [0.166666666666667; 0.166666666666667; 0.666666666666667]
	elseif(order==3)
		npoints = 4
		weights .= [-0.974278579257493; 0.902109795608790; 0.902109795608790; 0.902109795608790]
		coords1 .= [ 0.333333333333333; 0.600000000000000; 0.200000000000000; 0.200000000000000]
		coords2 .= [ 0.333333333333333; 0.200000000000000; 0.600000000000000; 0.200000000000000]
		coords3 .= [ 0.333333333333333; 0.200000000000000; 0.200000000000000; 0.600000000000000]
	elseif(order==4)
		npoints = 6
		weights .= [0.386908262797819; 0.386908262797819; 0.386908262797819; 0.190442006391807; 0.190442006391807; 0.190442006391807]
		coords1 .= [0.108103018168070; 0.445948490915965; 0.445948490915965; 0.816847572980459; 0.091576213509771; 0.091576213509771]
		coords2 .= [0.445948490915965; 0.108103018168070; 0.445948490915965; 0.091576213509771; 0.816847572980459; 0.091576213509771]
		coords3 .= [0.445948490915965; 0.445948490915965; 0.108103018168070; 0.091576213509771; 0.091576213509771; 0.816847572980459]
	else
		error("order ",order," not supported yet");
	end

	return nothing
end# }}}
