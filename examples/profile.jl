#!/Applications/Julia-1.6.app/Contents/Resources/julia/bin/julia 
using dJUICE
using MAT

#Load model from MATLAB file
file = matopen(joinpath(@__DIR__, "..", "data", "temp.mat"))
mat  = read(file, "md")
close(file)
md = model(mat)

md.timestepping.final_time = 2

#Solve transient simulation
@time md = solve(md, "tr")