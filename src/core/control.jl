using Enzyme
Enzyme.API.looseTypeAnalysis!(false)
Enzyme.API.strictAliasing!(false)
Enzyme.API.typeWarning!(false)

function Control_Core(md::model, femmodel::FemModel) #{{{
	# Compute gradient 
	computeGradient(md, femmodel)
end#}}}
function computeGradient(md::model, femmodel::FemModel) #{{{
	#independent variable
	α = md.inversion.independent
	#initialize derivative as 0
	∂J_∂α = zero(α)

	# zero ALL depth of the model, make sure we get correct gradient
	dfemmodel = Enzyme.Compiler.make_zero(Base.Core.Typeof(femmodel), IdDict(), femmodel)
	# compute the gradient
	println("CALLING AUTODIFF, prepare to die...")
	@time autodiff(Enzyme.Reverse, costfunction, Duplicated(femmodel, dfemmodel), Duplicated(α, ∂J_∂α))

	#Put gradient in results
	InputUpdateFromVectorx(femmodel, ∂J_∂α, GradientEnum, VertexSIdEnum)
	RequestedOutputsx(femmodel, [GradientEnum])
end#}}}
