using dJUICE
using Flux
using MAT
using BSON: @load
using StatsBase

file = matopen("/Users/gongcheng/Dartmouth/myJulia/DATA/Helheim_model.mat")
mat  = read(file, "md")
close(file)
md = model(mat)
md.friction.coefficient =  mat["friction"]["C"][:];
md.friction.p=3.0 .* ones(md.mesh.numberofvertices)
md.friction.q=zeros(md.mesh.numberofvertices)

#md=solve(md,"Stressbalance")


md2 = model(mat, useDNN=true)
@load "../data/Helheim_Weertman_NN.bson" nn dtx dty
md2.friction.dnnChain = nn;
md2.friction.coefficient =  mat["friction"]["C"][:];
md2.friction.dtx = dtx;
md2.friction.dty = dty;

md2=solve(md2,"Stressbalance")
