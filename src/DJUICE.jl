# Julia version of the Ice-sheet and Sea-level System Model
#
# Author: Mathieu Morlighem
# email:  mathieu.morlighem@dartmouth.edu

module DJUICE
const userdir = "./usr"
const coredir = "./core"

using Enzyme

function __init__()
    Enzyme.Compiler.RunAttributor[] = false
end

include("$userdir/classes.jl")
export model, WeertmanFriction, SchoofFriction, DNNFriction
include("$userdir/exp.jl")
export expread, ContourToNodes
include("$userdir/utils.jl")
export archread, issmdir
include("$userdir/triangle.jl")
export triangle
include("$userdir/triangle_issm.jl")
export triangle_issm
include("$userdir/parameterization.jl")
export setmask, InterpFromMeshToMesh2d, SetIceShelfBC
include("$coredir/solve.jl")
export solve
include("$userdir/plotmodel.jl")
export plotmodel
include("$userdir/dataprocess.jl")
export IntegrateOverDomain, GetAreas

end
