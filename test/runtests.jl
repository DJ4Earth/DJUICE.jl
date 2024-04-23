using DJUICE
using Test

function searchdir(path,key) 
	filter(x->occursin(key,x), readdir(path))
end

@time begin
	@time @testset "Model Struct Tests" begin include("modelstructtests.jl") end

	# test each individual cases, name with test[0-9]*.jl
	testsolutions = searchdir("./", r"test[0-9]*.jl")
	@time @testset "Model Solution Tests" begin
		for tf in testsolutions
			include(tf)
		end
	end

	# AD test
	@time include("testad.jl")
	#@time include("testad2.jl")

	# GPU test
	#@time include("testGPU.jl")
end
