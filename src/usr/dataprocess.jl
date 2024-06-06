function IntegrateOverDomain(md::model, data::Vector{Float64}) #{{{
	# IntegrateOverDomain - integrating data over the whole domain
	elements = md.mesh.elements
	x = md.mesh.x
	y = md.mesh.y

	# compute areas
	eleAreas=GetAreas(elements,x,y);

	eleData = 1/3*eleAreas.*(data[elements[:,1],:] + data[elements[:,2],:] + data[elements[:,3],:])
	return sum(eleData)
end#}}}
function GetAreas(index::Matrix{Int64}, x::Vector{Float64}, y::Vector{Float64}) #{{{
	# GetAreas - compute areas or volumes of elements
	nels=size(index,1);
	nods=length(x);
	areas=zeros(nels,1);

	x1=x[index[:,1]]; x2=x[index[:,2]]; x3=x[index[:,3]];
	y1=y[index[:,1]]; y2=y[index[:,2]]; y3=y[index[:,3]];

	return 0.5*((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1))
end#}}}
