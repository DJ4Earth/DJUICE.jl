#!/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia 
using dJUICE
using MAT
using Enzyme

include("./cost.jl")

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "..", "data","temp12k.mat"))
mat  = read(file, "md")
close(file)
md = model(mat)

# ========================================================

#define control
α = md.friction.coefficient

md.friction.coefficient = Active(α)

#initialize derivative as 0
∂J_∂α = zero(α)

#@show cost(md, α)

#Call enzyme to get derivative of cost function
Enzyme.API.looseTypeAnalysis!(true)
Enzyme.API.strictAliasing!(false)
autodiff(cost, Active, md, Duplicated(α, ∂J_∂α))
print(∂f_∂α[1:10])
