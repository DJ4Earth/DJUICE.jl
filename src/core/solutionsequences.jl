function solutionsequence_linear(femmodel::FemModel,analysis::Analysis) # {{{

	#First, update constraints in case the levelset has changed
	UpdateConstraintsx(femmodel, analysis)

	#Get new matrices
	Kff, Kfs, pf = SystemMatricesx(femmodel,analysis)

	#Enforce constraints
	ys = CreateNodalConstraintsx(femmodel.nodes)
	Reduceloadx!(pf, Kfs, ys)

	#Solve!
	uf = Solverx(Kff, pf)

	#Merge uf with ys and update inputs
	gsize = NumberOfDofs(femmodel.nodes,GsetEnum)
	ug = IssmVector(gsize)
	Mergesolutionfromftogx(ug, uf, ys, femmodel.nodes)

	InputUpdateFromSolutionx(analysis, ug, femmodel)

	return nothing
end# }}}
function solutionsequence_nonlinear(femmodel::FemModel,analysis::Analysis,maxiter::Int64,restol::Float64,reltol::Float64,abstol::Float64) # {{{

	#First, update constraints in case the levelset has changed
	UpdateConstraintsx(femmodel, analysis)

	#Initialize number of iterations
	count = 0
	converged = false

	#Get existing solution
	ug = GetSolutionFromInputsx(analysis,femmodel)
	uf = Reducevectorgtofx(ug,femmodel.nodes)

	#Update once again the solution to make sure that vx and vxold are similar (for next step in transient or steadystate)
	InputUpdateFromSolutionx(analysis,ug,femmodel)

	#Loop until we reach convergence
	while(~converged)

		#Get new matrices
		Kff, Kfs, pf = SystemMatricesx(femmodel,analysis)

		#Enforce constraints
		ys = CreateNodalConstraintsx(femmodel.nodes)
		Reduceloadx!(pf, Kfs, ys)

		#Solve!
		old_uf = uf
		uf = Solverx(Kff, pf, old_uf)

		#Merge uf with ys
		Mergesolutionfromftogx(ug, uf, ys, femmodel.nodes)

		#Check for convergence
		converged = convergence(Kff,pf,uf,old_uf,restol,reltol,abstol)
		InputUpdateFromSolutionx(analysis,ug,femmodel)

		#Increase count
		count += 1
		if(count>=maxiter)
			println("   maximum number of nonlinear iterations (",maxiter,") exceeded")
			converged = true
		end
	end

	print("\n   total number of iterations: ",  count,  "\n")
	return nothing

end# }}}
function convergence(Kff::IssmMatrix, pf::IssmVector, uf::IssmVector, old_uf::IssmVector, restol::Float64, reltol::Float64, abstol::Float64)#{{{

	print("   checking convergence\n");

	#If solution vector is empty, return true
	if(IsEmpty(uf))
		return true
	end

	#Convergence criterion #1: force equilibrium (Mandatory)
	#compute K[n]U[n-1] - F
	KUold  = Duplicate(uf);    MatMult!(Kff,old_uf,KUold)
	KUoldF = Duplicate(KUold); VecCopy!(KUold, KUoldF); AXPY!(KUoldF, -1.0, pf)
	nKUoldF = Norm(KUoldF,2)
	nF      = Norm(pf,2)
	res = nKUoldF/nF
	if ~isfinite(res)
		println("norm nf = ", nF, " and norm kuold = ",nKUoldF)
		error("mechanical equilibrium convergence criterion is not finite!")
	end
	if(res<restol)
		print("   mechanical equilibrium convergence criterion ", res*100, " < ", restol*100, " %\n")
		converged=true
	else
		print("   mechanical equilibrium convergence criterion ", res*100, " > ", restol*100, " %\n")
		converged=false;
	end

	#Convergence criterion #2: norm(du)/norm(u)
	if ~isnan(reltol)
		duf = Duplicate(old_uf); VecCopy!(old_uf,duf); AXPY!(duf, -1.0, uf)
		ndu = Norm(duf, 2); nu = Norm(old_uf, 2)
		if ~isfinite(ndu) | ~isfinite(nu) 
			error("convergence criterion is not finite!")
		end
		if((ndu/nu)<reltol)
			print("   Convergence criterion: norm(du)/norm(u)      ", ndu/nu*100, " < ", reltol*100, " %\n")
		else
			print("   Convergence criterion: norm(du)/norm(u)      ", ndu/nu*100, " > ", reltol*100, " %\n")
			converged=false;
		end
	end

	#Convergence criterion #3: max(du)
	if ~isnan(abstol)
		duf = Duplicate(old_uf); VecCopy!(old_uf,duf); AXPY!(duf, -1.0, uf)
		nduinf= Norm(duf, 3)
		if ~isfinite(nduinf) 
			error("convergence criterion is not finite!")
		end
		if(nduinf<abstol)
			print("   Convergence criterion: max(du)               ", nduinf, " < ", abstol, "\n")
		else
			print("   Convergence criterion: max(du)               ", nduinf, " > ", abstol, "\n")
			converged=false;
		end
	end

	return converged

end#}}}
