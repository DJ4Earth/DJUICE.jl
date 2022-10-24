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
md=solve(md, :Stressbalance)

mdnn = model(mat)
@load "../data/Helheim_friction_NN_bedonly.bson" nn dtx dty
mdnn = model(mdnn;friction=DNNFriction())

mdnn.friction.dnnChain = nn;
mdnn.friction.coefficient =  mat["friction"]["C"][:];
mdnn.friction.dtx = dtx;
mdnn.friction.dty = dty;
mdnn.geometry.ssx = mat["results"]["ssx"][:]
mdnn.geometry.ssy = mat["results"]["ssy"][:]
mdnn.geometry.bsx = mat["results"]["bsx"][:]
mdnn.geometry.bsy = mat["results"]["bsy"][:]

mdnn=solve(mdnn, :Stressbalance)

