#!/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia 
using dJUICE
using MAT
using Enzyme

include("./cost.jl")

#Load model from MATLAB file
file = matopen("./temp12k.mat")
mat  = read(file, "md")
close(file)
md = model(mat)

# ========================================================

#define control
α = md.friction.coefficient

md.friction.coefficient = Active(α)

#initialize derivative as 0
∂J_∂α = zero(α)

#Call enzyme to get derivative of cost function
autodiff(cost, Active, md, Duplicated(α, ∂J_∂α))
print(∂f_∂α[1:10])
