# Julia version of the Ice-sheet and Sea-level System Model
#
# Author: Mathieu Morlighem, Joel Wilner
# email:  mathieu.morlighem@dartmouth.edu

module dJUICE
const userdir = joinpath("..", "usr")
const coredir = joinpath("..", "core")

include("$userdir/classes.jl")
export model, WeertmanFriction
include("$userdir/exp.jl")
export expread, ContourToNodes
include("$userdir/utils.jl")
export archread, issmdir
include("$userdir/triangle.jl")
export triangle
include("$userdir/parameterization.jl")
export setmask, InterpFromMeshToMesh2d
include("$coredir/solve.jl")
export solve, solve2

end
