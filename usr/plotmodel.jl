import ColorSchemes.jet
using  GLMakie
using .ISSM

function plotmodel( md::model, data::Vector, showvertices::Bool=false, showfacets::Bool=true)

	vertexcolor  = :black
	facetcolor   = :blue

	if data isa AbstractVector

		if length(data)==md.mesh.numberofelements
			# vector of polygons
			x = md.mesh.x
			y = md.mesh.y
			index = md.mesh.elements
			ps = [Makie.GeometryBasics.Polygon([Point2(x[index[i,1]], y[index[i,1]]), Point2(x[index[i,2]], y[index[i,2]]), Point2(x[index[i,3]], y[index[i,3]])])
					for i in 1:md.mesh.numberofelements]

			fig, ax, h = Makie.poly(ps, color = data, colormap = jet)

			#Add colorbar
			Colorbar(fig[1, 2], limits = (minimum(data), maximum(data)), colormap = jet)
		elseif length(data)==md.mesh.numberofvertices
			fig, ax, h = Makie.mesh( [md.mesh.x md.mesh.y], md.mesh.elements, shading = false, color = data, colormap = jet)

			#Add colorbar
			#Colorbar(fig[1, 2], h, width=25)
		else
			error("data of size "*string(length(data))*" not supported yet!")
		end
	else
		# default to single color
		@assert length(data)==1
		fig, ax, h = Makie.mesh( [md.mesh.x md.mesh.y], md.mesh.elements, shading = false, color = data, colormap = jet)
	end

	if showfacets
		Makie.wireframe!(ax, h[1][], color=facetcolor)
	end

	if showvertices
		Makie.scatter!( [md.mesh.x md.mesh.y], markersize = 4, color = vertexcolor)
	end

	return fig
end

function plotmodel(md::model,data::BitVector) #{{{

	println("Converting BitVector to Vector")
	data2 = Vector{Float64}(undef,size(data))
	for i in 1:length(data)
		data2[i] = Float64(data[i])
	end

	plotmodel(md,data2)

end#}}}
function plotmodel(md::model,data::String) #{{{

	if(data=="mesh")
		poly([md.mesh.x md.mesh.y], md.mesh.elements, strokewidth=1, shading=false)
	else
		error(data, " plot not supported yet")
	end

end#}}}
