function Matrix2x2Determinant(A::Matrix{Float64}) #{{{

	return  A[1,1]*A[2,2]-A[2,1]*A[1,2]

end#}}}
function Matrix2x2Invert(A::Matrix{Float64}) #{{{

	#Initialize output
	Ainv = Matrix{Float64}(undef,2,2)

	#Compute determinant
	det = Matrix2x2Determinant(A)
	if(abs(det)<eps(Float64)) error("Determinant smaller than machine epsilon") end

	#Multiplication is faster than divsion, so we multiply by the reciprocal
	det_reciprocal = 1/det

	#compute invert matrix
	Ainv[1,1]=   A[2,2]*det_reciprocal # =  d/det
   Ainv[1,2]= - A[1,2]*det_reciprocal # = -b/det
   Ainv[2,1]= - A[2,1]*det_reciprocal # = -c/det
   Ainv[2,2]=   A[1,1]*det_reciprocal # =  a/det

	return  Ainv

end#}}}
