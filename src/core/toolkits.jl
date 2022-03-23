
#Matrix

#Toolkit #1: serial sparse arrays
#= DEACTIVATED FOR NOW
using SparseArrays
mutable struct IssmMatrix #{{{
	M::Int64
	N::Int64
	rows::Vector{Int64}
	cols::Vector{Int64}
	vals::Vector{Float64}
	matrix::SparseMatrixCSC{Float64,Int64}
end #}}}
function IssmMatrix(M::Int64,N::Int64)#{{{
	return IssmMatrix(M, N, Vector{Int64}(undef,0), Vector{Int64}(undef,0), Vector{Float64}(undef,0), spzeros(0,0))
end#}}}
function AddValues!(matrix::IssmMatrix,m::Int64,midx::Vector{Int64},n::Int64,nidx::Vector{Int64},values::Matrix{Float64})#{{{

	#This is inefficient now, but it will work
	for i in 1:m
		if(midx[i]==-1) continue end
		for j in 1:n
			if(nidx[j]==-1) continue end
			push!(matrix.rows, midx[i])
			push!(matrix.cols, nidx[j])
			push!(matrix.vals, values[i,j])
		end
	end

	return nothing
end#}}}
function GetSize(matrix::IssmMatrix)#{{{

	return size(matrix.matrix)

end#}}}
function Assemble!(matrix::IssmMatrix)#{{{

	matrix.matrix = sparse(matrix.rows, matrix.cols, matrix.vals, matrix.M, matrix.N)

	return nothing
end#}}}
=#

#Toolkit #2: dense matrix (for ensyme)
mutable struct IssmMatrix #{{{
	M::Int64
	N::Int64
	matrix::Matrix{Float64}
end #}}}
function IssmMatrix(M::Int64,N::Int64)#{{{
	return IssmMatrix(M, N, zeros(M,N))
end#}}}
function AddValues!(matrix::IssmMatrix,m::Int64,midx::Vector{Int64},n::Int64,nidx::Vector{Int64},values::Matrix{Float64})#{{{

	#This is inefficient now, but it will work
	for i in 1:m
		if(midx[i]==-1) continue end
		for j in 1:n
			if(nidx[j]==-1) continue end
			matrix.matrix[midx[i],nidx[j]] += values[i,j]
		end
	end

	return nothing
end#}}}
function GetSize(matrix::IssmMatrix)#{{{

	return size(matrix.matrix)

end#}}}
function Assemble!(matrix::IssmMatrix)#{{{

	#Nothing to do here :)
	return nothing
end#}}}


#Vector
mutable struct IssmVector #{{{
	vector::Vector{Float64}
end #}}}
function IssmVector(M::Int64)#{{{
	return IssmVector(zeros(M))
end#}}}
function GetSize(vector::IssmVector)#{{{

	return length(vector.vector)

end#}}}
function AddValues!(vector::IssmVector,m::Int64,midx::Vector{Int64},values::Vector{Float64})#{{{

	#This is inefficient now, but it will work
	for i in 1:m
		if(midx[i]==-1) continue end
		vector.vector[midx[i]] += values[i]
	end

	return nothing
end#}}}
function SetValues!(vector::IssmVector,m::Int64,midx::Vector{Int64},values::Vector{Float64})#{{{

	#This is inefficient now, but it will work
	for i in 1:m
		if(midx[i]==-1) continue end
		vector.vector[midx[i]] = values[i]
	end

	return nothing
end#}}}
function IsEmpty(vector::IssmVector)#{{{

	return GetSize(vector)==0

end#}}}
function Duplicate(vector::IssmVector)#{{{

	#Copy data structure
	M=GetSize(vector)
	return IssmVector(M)

end#}}}
function VecCopy!(x::IssmVector,y::IssmVector)#{{{

	y.vector = x.vector

	return nothing
end#}}}
function Assemble!(vector::IssmVector)#{{{

	#Nothing to do for this toolkit
	return nothing

end#}}}
function ToSerial(vector::IssmVector)#{{{

	return vector.vector

end#}}}
function Norm(x::IssmVector,type::Int64)#{{{

	norm = 0.0

	if type==2
		for i in 1:length(x.vector)
			norm += x.vector[i]^2
		end
		norm = sqrt(norm)
	elseif type==3
		#Infinite norm
		for i in 1:length(x.vector)
			if(abs(x.vector[i])>norm) norm = abs(x.vector[i]) end
		end
	else
		error("type ",type," not supported yet")
	end

	return norm

end#}}}

#Operations
function MatMult!(A::IssmMatrix,x::IssmVector,y::IssmVector) #{{{

	y.vector = A.matrix*x.vector

	return nothing
end#}}}
function AXPY!(y::IssmVector,alpha::Float64,x::IssmVector) #{{{

	y.vector = alpha*x.vector + y.vector

	return nothing
end#}}}
function Solverx(A::IssmMatrix, b::IssmVector, xold::IssmVector) #{{{

	#Initialize output
	#x = IssmVector(GetSize(xold))
	
	return Solverx(A, b)

end#}}}
function Solverx(A::IssmMatrix, b::IssmVector) #{{{

	#Initialize output
	x = IssmVector(0)

	#Solve linear system
	x.vector = A.matrix\b.vector

	return x

end#}}}
