
#define function
function cost(md, α)

	#change md according to incoming α
	md.friction.coefficient = α
	md.inversion.iscontrol  = true

	#Get cost function
	J = solve2(md)

	return J
end
