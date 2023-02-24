using dJUICE
using Flux
using MAT
using BSON: @load
using StatsBase
using Random

# file all files in the folder

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
mdnn.geometry.ssx = mat["results"]["ssx"][:]
mdnn.geometry.ssy = mat["results"]["ssy"][:]
mdnn.geometry.bsx = mat["results"]["bsx"][:]
mdnn.geometry.bsy = mat["results"]["bsy"][:]
mdnn.stressbalance.restol = 0.02;
mdnn.stressbalance.reltol = 0.2;

mdnn = model(mdnn;friction=DNNFriction())
# load all the DNNs in the folder
fileList =  readdir(pwd()*"/../DATA/meanDNNs/")

for i in 1:length(fileList)
	@load pwd()*"/../DATA/Helheim_friction_NN_Schoof.bson" nn dtx dty
	push!(mdnn.friction.dnnChain, nn)
	push!(mdnn.friction.dtx, dtx)
	push!(mdnn.friction.dty, dty)
end

mdnn=solve(mdnn, :Stressbalance)

