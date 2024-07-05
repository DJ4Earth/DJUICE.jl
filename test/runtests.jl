using DJUICE
using Test

function searchdir(path,key) 
	filter(x->occursin(key,x), readdir(path))
end

function compareArchive(tf::String, procedure::Symbol)
	# run test
	@inline include(tf)
	id = match(r"test(\d+).jl", tf).captures[1]
	archive_name = "Archive"*string(id)
	archive_path = issmdir()*"/test/Archives/"*archive_name*".arch"
	if procedure===:update
		# update Archive
	else
		# check Archive
		if isfile(archive_path)
			for k=1:length(field_names)
				# Get field and tolerance
				field=field_values[k];
				fieldname=field_names[k];
				tolerance=field_tolerances[k];

				# Compare to archive
				if !isnan(tolerance)
					# Our output is in the correct order (n,1) or (1,1), so we do not need to transpose again
					archive = archread(archive_path, archive_name*"_field"*string(k))
					error_diff = (maximum(abs.(archive-field))/(maximum(abs.(archive))+eps(Float64)))

					@test isnan(error_diff) == false
					@test error_diff < tolerance skip = isnan(tolerance)
				end
			end
		else
			@warn "$archive_name does not exist! Skip the comparison of the results"
		end
	end
end

@time @testset "DJUICE" begin
	@time @testset "Model Struct Tests" begin include("modelstructtests.jl") end

	# test each individual cases, name with test[0-9]*.jl
	testsolutions = searchdir("./", r"test\d+.jl")
	for tf in testsolutions
		@time @testset "Model Solution Tests: $tf"  begin
			# check the results vs. saved archive
			compareArchive(tf, :test)
		end
	end

	# GPU test
	#@time include("testGPU.jl")
end
