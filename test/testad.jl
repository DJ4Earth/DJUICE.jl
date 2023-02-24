#!/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia --project
using dJUICE
using MAT
using Enzyme
using Enzyme_jll

#define cost function
function cost(md::model, frictioncoeff::Vector{Float64})
    
    #Set friction coefficient based on input
    md.friction.coefficient = frictioncoeff
    
    #Solve stress balance
    md = solve(md, :Stressbalance)
    
    #return misfit to observations
    vel_data  = sqrt.(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2)
    vel_model = md.results["StressbalanceSolution"]["Vel"]
    return sum(sqrt.((vel_data - vel_model).^2))*md.constants.yts
end

#Load model from MATLAB file
#file = matopen(joinpath(@__DIR__, "..", "data","temp12k.mat")) #BIG model
file = matopen(joinpath(@__DIR__, "..", "data","temp.mat")) #SMALL model (35 elements)
mat  = read(file, "md")
close(file)
md = model(mat)

#make model run faster 
md.stressbalance.maxiter = 1

#Call cost once to compile it
#@time println("\n\nInitial cost function is J = ", cost(md, md.friction.coefficient), " m/yr (1st call)")

#Call cost again to test how long it takes to run
#@time println("\n\nInitial cost function is J = ", cost(md, md.friction.coefficient), " m/yr (2d call)")

#Now call AD!
md.inversion.iscontrol = 1
md = solve2(md, true)
#print(∂f_∂α[1:10])
