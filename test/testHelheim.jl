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
@load "../data/Helheim_friction_NN_Schoof.bson" nn dtx dty
mdnn = model(mdnn;friction=DNNFriction())

mdnn.friction.dnnChain = nn;
mdnn.friction.coefficient =  mat["friction"]["C"][:];
mdnn.friction.dtx = dtx;
mdnn.friction.dty = dty;
mdnn.friction.Cmax = 0.8
mdnn.friction.velThreshold = 10e6./mdnn.constants.yts
mdnn.geometry.ssx = mat["results"]["ssx"][:]
mdnn.geometry.ssy = mat["results"]["ssy"][:]
mdnn.geometry.bsx = mat["results"]["bsx"][:]
mdnn.geometry.bsy = mat["results"]["bsy"][:]
mdnn.stressbalance.restol = 0.01;
mdnn.stressbalance.reltol = NaN;

mdnn=solve(mdnn, :Stressbalance)

