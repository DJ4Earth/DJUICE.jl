#!/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia --project
using dJUICE
using MAT
using Enzyme

include("./cost.jl")

#Load model from MATLAB file
#file = matopen(joinpath(@__DIR__, "..", "data","temp12k.mat")) #BIG model
file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMAL model (35 elements)
mat  = read(file, "md")
close(file)

# ========================================================

md = model(mat)

#make model run faster 
md.stressbalance.maxiter = 1

# d_md = copy(md)

#define control
α = md.friction.coefficient
#initialize derivative as 0
∂J_∂α = zero(α)

#Call enzyme to get derivative of cost function
Enzyme.API.looseTypeAnalysis!(true)
Enzyme.API.strictAliasing!(false)
# TODO: We might have to make this `Duplicated(md, d_md)`
# TODO(@wsmoses): How do we make this sparsely active?
#                 We could construct the model as part of the code to differentiate...
autodiff(cost, Active, md, Duplicated(α, ∂J_∂α))
print(∂f_∂α[1:10])
