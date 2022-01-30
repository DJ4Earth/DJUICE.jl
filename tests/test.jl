#!/Applications/Julia-1.6.app/Contents/Resources/julia/bin/julia 
include("../issm.jl")
using .ISSM
using MAT

#Load model from MATLAB file
#file = matopen("./temp.mat")
#file = matopen("/Users/mmorligh/Desktop/issmuci/trunk-jpl/test/NightlyRun/temp.mat")
#file = matopen("/Users/mmorligh/Desktop/issmuci/trunk-jpl/examples/PigSensitivity/temp.mat")
#file = matopen("/Users/mmorligh/Desktop/issmuci/trunk-jpl/examples/PigSensitivity/temp12k.mat")
mat  = read(file, "md")
close(file)
md = model(mat)

println("NUMBEROFELEMENTS: ",md.mesh.numberofelements)
md.timestepping.final_time = 2

#Solve stress balance
@time md = solve(md, "tr")
