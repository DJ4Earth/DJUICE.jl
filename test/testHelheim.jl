using dJUICE
using Flux
using MAT
using BSON: @load

file = matopen("/Users/gongcheng/Dartmouth/myJulia/DATA/Helheim_model.mat")
mat  = read(file, "md")
close(file)
md = model(mat, useDNN=true)
@load "../data/Helheim_Weertman_NN.bson" nn C₀ dtx dty
md.friction.dnnChain = nn;
md.friction.coefficient =  mat["friction"]["C"][:];
md.friction.dtx = dtx;
md.friction.dty = dty;

md=solve(md,"Stressbalance")
