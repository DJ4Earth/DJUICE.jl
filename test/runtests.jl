using DJUICE
using Test

include("utils.jl")

@time @testset "DJUICE" begin
	@time @testset "Model Struct Tests" begin include("modelstructtests.jl") end

	# test each individual cases, name with test[0-9]*.jl
	testsolutions = searchdir(joinpath(dirname(pathof(DJUICE)), "../test/"), r"test\d+.jl")
	for tf in testsolutions
		@time @testset "Model Solution Tests: $tf"  begin
			# check the results vs. saved archive
			include(tf)
		end
	end
end
