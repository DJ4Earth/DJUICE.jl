using dJUICE
using Flux
using MAT
using BSON: @load
using StatsBase
using Random

file = matopen(pwd()*"/../DATA/Helheim_model.mat")
mat  = read(file, "md")
close(file)
md = model(mat)
md = model(md;friction=dJUICE.SchoofFriction())
md.friction.C =  mat["friction"]["C"][:];
md.friction.m =  mat["friction"]["m"][:];
md.friction.Cmax =  mat["friction"]["Cmax"][:];
md=solve(md, :Stressbalance)

mdnn = model(mat)
@load pwd()*"/../DATA/Helheim_friction_NN_Schoof.bson" nn dtx dty
mdnn = model(mdnn;friction=DNNFriction())

mdnn.friction.dnnChain = nn;
mdnn.friction.coefficient =  mat["friction"]["C"][:];
mdnn.friction.dtx = dtx;
mdnn.friction.dty = dty;
mdnn.friction.Cmax = 0.5
mdnn.friction.velThreshold = 100e6./mdnn.constants.yts
mdnn.geometry.ssx = mat["results"]["ssx"][:]
mdnn.geometry.ssy = mat["results"]["ssy"][:]
mdnn.geometry.bsx = mat["results"]["bsx"][:]
mdnn.geometry.bsy = mat["results"]["bsy"][:]
mdnn.stressbalance.restol = 0.02;
mdnn.stressbalance.reltol = 0.2;

mdnn=solve(mdnn, :Stressbalance)

