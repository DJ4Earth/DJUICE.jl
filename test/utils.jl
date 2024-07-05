using DJUICE
using Test

function searchdir(path,key)
	filter(x->occursin(key,x), readdir(path))
end

function compareArchive(tf::String, field_names::Vector{String}, field_tolerances::Vector{Float64}, field_values::Vector{Vector{Float64}}, procedure::Symbol)
	# make the test
	id = match(r"test(\d+).jl", tf).captures[1]
	archive_name = "Archive"*string(id)
	archive_path = issmdir()*"/test/Archives/"*archive_name*".arch"
	if procedure===:update
		# TODO: update Archive
	else
		# check Archive
		if isfile(archive_path)
			@time @testset "    Compare with $archive_name"  begin
				for k=1:length(field_names)
					# Compare to archive
					if !isnan(field_tolerances[k])
						@time @testset "      $(field_names[k]): "  begin
							# Our output is in the correct order (n,1) or (1,1), so we do not need to transpose again
							archive = archread(archive_path, archive_name*"_field"*string(k))
							error_diff = (maximum(abs.(archive-field_values[k]))/(maximum(abs.(archive))+eps(Float64)))

							@test isnan(error_diff) == false
							@test error_diff < field_tolerances[k] skip = isnan(field_tolerances[k])
						end
					end
				end
			end
		else
			@warn "$archive_name does not exist! Skip the comparison of the results"
		end
	end
end
