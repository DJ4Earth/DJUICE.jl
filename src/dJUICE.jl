# Julia version of the Ice-sheet and Sea-level System Model
#
# Author: Mathieu Morlighem
# email:  mathieu.morlighem@dartmouth.edu

module dJUICE
const userdir = "./usr"
const coredir = "./core"

include("$userdir/classes.jl")
export model, WeertmanFriction
include("$userdir/exp.jl")
export expread, ContourToNodes
include("$userdir/utils.jl")
export archread, issmdir
include("$userdir/triangle.jl")
export triangle
include("$userdir/triangle2.jl")
export triangle2
include("$userdir/parameterization.jl")
export setmask, InterpFromMeshToMesh2d
include("$coredir/solve.jl")
export solve, solve2

end
