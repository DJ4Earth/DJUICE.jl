#Matice class definition
struct Matice#{{{
	vx_input::Input
	vy_input::Input
	B_input::Input
	n_input::Input
end# }}}

function Matice(element::Tria) #{{{

	vx_input  = GetInput(element, VxEnum)
	vy_input  = GetInput(element, VyEnum)
	B_input   = GetInput(element, MaterialsRheologyBEnum)
	n_input   = GetInput(element, MaterialsRheologyNEnum)

	return Matice(vx_input, vy_input, B_input, n_input)
end#}}}

#vertices functions
function ViscositySSA(matice::Matice, xyz_list::Matrix{Float64}, gauss::GaussTria, i::Int64) #{{{

	#Get strain rate
	dvx = GetInputDerivativeValue(matice.vx_input,xyz_list,gauss,i)
	dvy = GetInputDerivativeValue(matice.vy_input,xyz_list,gauss,i)
	eps_xx = dvx[1]
	eps_yy = dvy[2]
	eps_xy = 0.5*(dvx[2] + dvy[1])

	#In SSA, eps_eff^2 = exx^2 + eyy^2 + exy^2 + exx*eyy
	eps_eff = sqrt(eps_xx*eps_xx + eps_yy*eps_yy + eps_xy*eps_xy + eps_xx*eps_yy)

	#Get B and n
	n = GetInputValue(matice.n_input, gauss, i)
	B = GetInputValue(matice.B_input, gauss, i)

	#Compute viscosity
	if eps_eff==0.
		mu = 1.e+14/2
	else
		mu = B/(2*eps_eff^((n-1)/n))
	end

	return mu::Float64
end #}}}
