#Node class definition
mutable struct Node #{{{
	id::Int64
	sid::Int64
	indexingupdate::Bool
	active::Bool
	gsize::Int64
	gdoflist::Vector{Int64}
	fdoflist::Vector{Int64}
	sdoflist::Vector{Int64}
	svalues::Vector{Float64}
end# }}}

#Node functions
function Base.show(io::IO, this::Node)# {{{

	println(io,"Node:")
	println(io,"   id:  ",this.id)
	println(io,"   sid: ",this.sid)
	println(io,"   indexingupdate: ",this.indexingupdate)
	println(io,"   gsize: ",this.gsize)
	println(io,"   gdoflist: ",this.gdoflist)
	println(io,"   fdoflist: ",this.fdoflist)
	println(io,"   sdoflist: ",this.sdoflist)
	println(io,"   svalues: ",this.svalues)
end# }}}
function Activate!(node::Node) #{{{

	if(!node.active)
		node.indexingupdate = true
		node.active = true
		for i in 1:node.gsize
			node.fdoflist[i]  = +1
			node.sdoflist[i]  = -1
			node.svalues[i]   = 0.0
		end
	end

end# }}}
function Deactivate!(node::Node) #{{{

	if(node.active)
		node.indexingupdate = true
		node.active = false
		for i in 1:node.gsize
			node.fdoflist[i]  = -1
			node.sdoflist[i]  = +1
			node.svalues[i]   = 0.0
		end
	end

end# }}}
function ApplyConstraint(node::Node,dof::Int8,value::Float64) #{{{

	node.indexingupdate = true
	node.fdoflist[dof]  = -1
	node.sdoflist[dof]  = +1
	node.svalues[dof]   = value

end# }}}
function CreateNodalConstraints(node::Node,ys::IssmVector) #{{{

	if(SSize(node)>0)
		SetValues!(ys,node.gsize,node.sdoflist,node.svalues)
	end

end# }}}
function DistributeDofs(node::Node,setenum::IssmEnum,dofcount::Int64) #{{{

	if setenum==GsetEnum
		for i in 1:node.gsize
			node.gdoflist[i] = dofcount
			dofcount += 1
		end
	elseif setenum==FsetEnum
		for i in 1:node.gsize
			if  node.fdoflist[i]!=-1
				@assert node.sdoflist[i]==-1
				node.fdoflist[i] = dofcount
				dofcount += 1
			end
		end
	elseif setenum==SsetEnum
		for i in 1:node.gsize
			if  node.sdoflist[i]!=-1
				@assert node.fdoflist[i]==-1
				node.sdoflist[i] = dofcount
				dofcount += 1
			end
		end
	else
		error("not supported")
	end

	return dofcount
end# }}}
function GetNumberOfDofs(node::Node,setenum::IssmEnum) #{{{

	if setenum==GsetEnum
		dofcount = node.gsize
	elseif setenum==FsetEnum
		dofcount = 0
		for i=1:node.gsize
			if  node.fdoflist[i]!=-1
				dofcount += 1
			end
		end
	elseif setenum==SsetEnum
		dofcount = 0
		for i=1:node.gsize
			if  node.sdoflist[i]!=-1
				dofcount += 1
			end
		end
	else
		error("not supported")
	end

	return dofcount

end# }}}
function GetDofList(node::Node,doflist::Vector{Int64},count::Int64,setenum::IssmEnum) #{{{

	if setenum==GsetEnum
		for i in 1:node.gsize
			count += 1
			doflist[count] = node.gdoflist[i]
		end
	elseif setenum==FsetEnum
		for i=1:node.gsize
			#if  node.fdoflist[i]!=-1
				count += 1
				doflist[count] = node.fdoflist[i]
			#end
		end
	elseif setenum==SsetEnum
		for i=1:node.gsize
			#if  node.sdoflist[i]!=-1
				count += 1
				doflist[count] = node.sdoflist[i]
			#end
		end
	else
		error("not supported")
	end

	return count

end# }}}
function GetGlobalDofList(nodes::Vector{Node},ndofs::Int64,setenum::IssmEnum) #{{{

	#Allocate list
	doflist = Vector{Int64}(undef,ndofs)

	#Assign values 
	count = 0
	for i in 1:length(nodes)
		count = GetDofList(nodes[i],doflist,count,setenum)
	end
	@assert count==ndofs

	return doflist

end# }}}
function VecReduce(node::Node,ug::Vector{Float64},uf::IssmVector) #{{{

	for i=1:node.gsize
		if node.fdoflist[i]!=-1
			uf.vector[node.fdoflist[i]] = ug[node.gdoflist[i]]
		end
	end

end# }}}
function VecMerge(node::Node,ug::IssmVector,uf::Vector{Float64},ys::Vector{Float64}) #{{{

	fsize = FSize(node)
	ssize = SSize(node)

	if fsize>0
		indices = Vector{Int64}(undef,fsize)
		values  = Vector{Float64}(undef,fsize)

		count = 1
		for i=1:node.gsize
			if node.fdoflist[i]!=-1
				indices[count] = node.gdoflist[i]
				values[count]  = uf[node.fdoflist[i]]
				count += 1
			end
		end
		SetValues!(ug,fsize,indices,values)
	end

	if ssize>0
		indices = Vector{Int64}(undef,ssize)
		values  = Vector{Float64}(undef,ssize)

		count = 1
		for i=1:node.gsize
			if node.sdoflist[i]!=-1
				indices[count] = node.gdoflist[i]
				values[count]  = ys[node.sdoflist[i]]
				count += 1
			end
		end
		SetValues!(ug,ssize,indices,values)
	end

end# }}}
function SSize(node::Node) #{{{

	ssize = 0

	for i=1:node.gsize
		if node.sdoflist[i]!=-1
			ssize+=1
		end
	end

	return ssize

end# }}}
function FSize(node::Node) #{{{

	fsize = 0

	for i=1:node.gsize
		if node.fdoflist[i]!=-1
			fsize+=1
		end
	end

	return fsize

end# }}}

#Nodes functions
function RequiresDofReindexing(nodes::Vector{Node}) #{{{

	for i in 1:length(nodes)
		if nodes[i].indexingupdate
			return true
		end
	end

	return false

end# }}}
function DistributeDofs(nodes::Vector{Node},setenum::IssmEnum) #{{{

	dofcount = 1

	for i in 1:length(nodes)
		dofcount = DistributeDofs(nodes[i],setenum,dofcount)
	end


end# }}}
function NumberOfDofs(nodes::Vector{Node},setenum::IssmEnum) #{{{

	numdofs = 0
	for i in 1:length(nodes)
		numdofs += GetNumberOfDofs(nodes[i],setenum)
	end
	return numdofs

end# }}}
