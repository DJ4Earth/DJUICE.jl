using dJUICE
using Flux
using MAT
using BSON: @load
using StatsBase
using Random

file = matopen("/Users/gongcheng/Dartmouth/myJulia/DATA/Helheim_model.mat")
mat  = read(file, "md")
close(file)
md = model(mat)
md = model(md;friction=dJUICE.SchoofFriction())
md.friction.C =  mat["friction"]["C"][:];
md.friction.m =  mat["friction"]["m"][:];
md.friction.Cmax =  mat["friction"]["Cmax"][:];
md=solve(md, :Stressbalance)

mdnn = model(mat)
@load "../data/Helheim_friction_NN_Schoof.bson" nn dtx dty
mdnn = model(mdnn;friction=DNNFriction())

mdnn.friction.dnnChain = nn;
mdnn.friction.coefficient =  mat["friction"]["C"][:];
mdnn.friction.dtx = dtx;
mdnn.friction.dty = dty;
mdnn.friction.Cmax = 0.8
mdnn.friction.velThreshold = 30e5./mdnn.constants.yts
mdnn.geometry.ssx = mat["results"]["ssx"][:]
mdnn.geometry.ssy = mat["results"]["ssy"][:]
mdnn.geometry.bsx = mat["results"]["bsx"][:]
mdnn.geometry.bsy = mat["results"]["bsy"][:]
mdnn.stressbalance.restol = 0.01;
mdnn.stressbalance.reltol = NaN;

mdnn=solve(mdnn, :Stressbalance)

