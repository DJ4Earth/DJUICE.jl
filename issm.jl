# Julia version of the Ice-sheet and Sea-level System Model
#
# Author: Mathieu Morlighem, Joel Wilner
# email:  mathieu.morlighem@dartmouth.edu

module ISSM

include("usr/classes.jl")
export model, WeertmanFriction
include("usr/exp.jl")
export expread, ContourToNodes
include("usr/utils.jl")
export archread, issmdir
include("usr/triangle.jl")
export triangle
include("usr/parameterization.jl")
export setmask, InterpFromMeshToMesh2d
include("core/solve.jl")
export solve, solve2

end
