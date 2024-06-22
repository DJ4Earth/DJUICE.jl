using DJUICE
using Test

function searchdir(path,key) 
	filter(x->occursin(key,x), readdir(path))
end

function compareArchive(id, procedure::Symbol)
	archive_name = "Archive"*string(id)
	archive_path = issmdir()*"/test/Archives/"*archive_name*".arch"
	if procedure===:update
		# update Archive
	else
		# check Archive
		if isempty(archive_path)
			@warn "$archive_name does not exist! Skip the comparison of the results"
		else
			for k=1:length(field_names)
				# Get field and tolerance
				field=field_values[k];
				fieldname=field_names[k];
				tolerance=field_tolerances[k];

				# Compare to archive
				# Our output is in the correct order (n,1) or (1,1), so we do not need to transpose again
				archive = archread(archive_path, archive_name*"_field"*string(k))
				error_diff = (maximum(abs.(archive-field))/(maximum(abs.(archive))+eps(Float64)))

				@test isnan(error_diff) == false
				@test error_diff < tolerance
			end
		end
	end
end

@time begin
	@time @testset "Model Struct Tests" begin include("modelstructtests.jl") end

	# test each individual cases, name with test[0-9]*.jl
	testsolutions = searchdir("./", r"test\d+.jl")
	@time @testset "Model Solution Tests" begin
		for tf in testsolutions
			# run the test
			include(tf)
			# check the results vs. saved archive
			testid = match(r"test(\d+).jl", tf).captures[1]
			compareArchive(testid, :test)
		end
	end

	# GPU test
	#@time include("testGPU.jl")
end
