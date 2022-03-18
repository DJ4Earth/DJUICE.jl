mutable struct ElementMatrix#{{{
	nrows::Int64
	gglobaldoflist::Vector{Int64}
	fglobaldoflist::Vector{Int64}
	sglobaldoflist::Vector{Int64}
	values::Matrix{Float64}
end #}}}
function ElementMatrix(nodes::Vector{Node})#{{{

	#Get matrix size
	nrows = NumberOfDofs(nodes,GsetEnum)

	#Initialize element matrix with zeros
	values = zeros(nrows,nrows)

	#Get dof lists
	gglobaldoflist=GetGlobalDofList(nodes,nrows,GsetEnum)
	fglobaldoflist=GetGlobalDofList(nodes,nrows,FsetEnum)
	sglobaldoflist=GetGlobalDofList(nodes,nrows,SsetEnum)

	return ElementMatrix(nrows,gglobaldoflist,fglobaldoflist,sglobaldoflist,values)
end#}}}
function Base.show(io::IO, this::ElementMatrix)# {{{

	println(io,"ElementMatrix:")
	println(io,"   nrows: ",this.nrows)
	println(io,"   gglobaldoflist: ",this.gglobaldoflist)
	println(io,"   fglobaldoflist: ",this.fglobaldoflist)
	println(io,"   sglobaldoflist: ",this.sglobaldoflist)
	print(io,"   values: ")
	display(this.values)

	return nothing
end# }}}
function AddToGlobal!(Ke::ElementMatrix,Kff::IssmMatrix,Kfs::IssmMatrix)#{{{

	#First check that the element matrix looks alright
	CheckConsistency(Ke)

	#See if we need to do anything
	is_fset = false
	is_sset = false
	for i in 1:Ke.nrows
		if(Ke.fglobaldoflist[i]>0) is_fset = true end
		if(Ke.sglobaldoflist[i]>0) is_sset = true end
	end

	if is_fset
		AddValues!(Kff,Ke.nrows,Ke.fglobaldoflist,Ke.nrows,Ke.fglobaldoflist,Ke.values)
	end
	if is_sset
		AddValues!(Kfs,Ke.nrows,Ke.fglobaldoflist,Ke.nrows,Ke.sglobaldoflist,Ke.values)
	end

	return nothing
end#}}}
function CheckConsistency(Ke::ElementMatrix)#{{{

	for i in 1:Ke.nrows
		for j in 1:Ke.nrows
			if(isnan(Ke.values[i,j])) error("NaN found in Element Matrix") end
			if(isinf(Ke.values[i,j])) error("Inf found in Element Matrix") end
			if(abs(Ke.values[i,j])>1.e+50) error("Element Matrix values exceeds 1.e+50") end
		end
	end

	return nothing
end#}}}

mutable struct ElementVector#{{{
	nrows::Int64
	fglobaldoflist::Vector{Int64}
	values::Vector{Float64}
end #}}}
function ElementVector(nodes::Vector{Node})#{{{

	#Get matrix size
	nrows = NumberOfDofs(nodes,GsetEnum)

	#Initialize element matrix with zeros
	values = zeros(nrows)

	#Get dof list
	fglobaldoflist=GetGlobalDofList(nodes,nrows,FsetEnum)

	return ElementVector(nrows,fglobaldoflist,values)
end#}}}
function Base.show(io::IO, this::ElementVector)# {{{

	println(io,"ElementVector:")
	println(io,"   nrows: ",this.nrows)
	println(io,"   fglobaldoflist: ",this.fglobaldoflist)
	print(io,"   values: ")
	display(this.values)

	return nothing
end# }}}
function AddToGlobal!(pe::ElementVector,pf::IssmVector)#{{{

	#First check that the element matrix looks alright
	CheckConsistency(pe)

	#See if we need to do anything
	is_fset = false
	for i in 1:pe.nrows
		if(pe.fglobaldoflist[i]>0)
			is_fset = true 
			break
		end
	end

	if is_fset
		AddValues!(pf,pe.nrows,pe.fglobaldoflist,pe.values)
	end

	return nothing
end#}}}
function CheckConsistency(pe::ElementVector)#{{{

	for i in 1:pe.nrows
		if(isnan(pe.values[i])) error("NaN found in Element Vector") end
		if(isinf(pe.values[i])) error("Inf found in Element Vector") end
		if(abs(pe.values[i])>1.e+50) error("Element Vector values exceeds 1.e+50") end
	end

	return nothing
end#}}}
