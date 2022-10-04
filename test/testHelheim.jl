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


md2 = model(mat)
@load "../data/Helheim_Weertman_NN.bson" nn dtx dty
md2 = model(md2;friction=DNNFriction())

md2.friction.dnnChain = nn;
md2.friction.coefficient =  mat["friction"]["C"][:];
md2.friction.dtx = dtx;
md2.friction.dty = dty;

#md2=solve(md2,"Stressbalance")

md3 = model(mat)
@load "../data/Helheim_friction_NN.bson" nn dtx dty
md3 = model(md3;friction=DNNFriction())

md3.friction.dnnChain = nn;
md3.friction.coefficient =  mat["friction"]["C"][:];
md3.friction.dtx = dtx;
md3.friction.dty = dty;

md3=solve(md3,"Stressbalance")
