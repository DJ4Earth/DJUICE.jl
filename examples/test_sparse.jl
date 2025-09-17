using Enzyme
Enzyme.API.strictAliasing!(false)
Enzyme.Compiler.VERBOSE_ERRORS[] = true
using DJUICE

n = 5
A = DJUICE.IssmMatrix(n,n)

for i = 1:n
	DJUICE.AddValues!(A, 1, [i], 1,[i], ones(1,1))
end
DJUICE.Assemble!(A)
b = DJUICE.IssmVector(n)
x = DJUICE.IssmVector(n)


function mybackslash(x, A, b) 
	x.vector = A.matrix \ b.vector
	norm(x.vector)
end


∂z_∂A = DJUICE.IssmMatrix(n,n)
∂z_∂b = DJUICE.IssmVector(n)
∂z_∂x = DJUICE.IssmVector(n)

Enzyme.autodiff(set_runtime_activity(Enzyme.Reverse), mybackslash,Active, Duplicated(x, ∂z_∂x), Duplicated(A, ∂z_∂A), Duplicated(b, ∂z_∂b))
