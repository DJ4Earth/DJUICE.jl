# DJUICE.jl
[![Build Status](https://github.com/DJ4Earth/DJUICE.jl/workflows/CI/badge.svg)](https://github.com/DJ4Earth/DJUICE.jl/actions)

**Differentiable JUlia ICE model**

## Overview

DJUICE.jl is a Julia package designed for differentiable ice sheet modeling, leveraging the power of [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) for automatic differentiation. This package employs the Finite Element method (FEM) and follows the structure of the [ice-sheet and Sea-level System Model (ISSM)](https://issm.jpl.nasa.gov/), which is the C++ counterpart for ice sheet modeling.

## Features

- Differentiable ice sheet modeling using Enzyme.jl.
- Finite Element method implementation for ice numerical sheet modeling.
- Based on the well-established ISSM framework.

## Installation

To install DJUICE.jl, you can use the Julia package manager. In your Julia REPL, run:

```julia
using Pkg
Pkg.add("DJUICE")
```

## Usage

Here is a simple example of how to use DJUICE.jl:

```julia
using DJUICE

# Initialize model
md = model()

# Create mesh
md = triangle(md,issmdir()*"/test/Exp/Square.exp",50000.)
md = setmask(md,"all","")

#Geometry
hmin=300.
hmax=1000.
ymin=minimum(md.mesh.y)
ymax=maximum(md.mesh.y)
xmin=minimum(md.mesh.x)
xmax=maximum(md.mesh.x)
md.geometry.thickness = hmax .+ (hmin-hmax)*(md.mesh.y .- ymin)./(ymax-ymin) .+ 0.1*(hmin-hmax)*(md.mesh.x .- xmin)./(xmax-xmin)
md.geometry.base      = -md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness
md.geometry.surface   = md.geometry.base+md.geometry.thickness
md.geometry.bed       = md.geometry.base .-10

#Initial velocity
md.initialization.vx=zeros(md.mesh.numberofvertices)
md.initialization.vy=zeros(md.mesh.numberofvertices)

# Physical parameters
md.materials.rheology_B=1.815730284801701e+08*ones(md.mesh.numberofvertices)
md.materials.rheology_n=3*ones(md.mesh.numberofelements)
md.friction.coefficient=20*ones(md.mesh.numberofvertices)
md.friction.p=ones(md.mesh.numberofvertices)
md.friction.q=ones(md.mesh.numberofvertices)

# Numerical tolerances
md.stressbalance.restol=0.05
md.stressbalance.reltol=0.05
md.stressbalance.abstol=NaN

#Boundary conditions
md.stressbalance.spcvx = NaN*ones(md.mesh.numberofvertices)
md.stressbalance.spcvy = NaN*ones(md.mesh.numberofvertices)
pos = findall(md.mesh.vertexonboundary)
md.stressbalance.spcvx[pos] .= 0.0
md.stressbalance.spcvy[pos] .= 0.0

# Solve
md=solve(md,:Stressbalance)
```

For detailed tutorials and examples, please refer to the [documentation](https://github.com/DJ4Earth/DJUICE.jl).

## Documentation

Comprehensive documentation is available to help you get started and make the most out of DJUICE.jl. It includes tutorials, API references, and examples. Access the documentation [here](https://github.com/DJ4Earth/DJUICE.jl).

## Contributing

Contributions are welcome! If you find a bug or have a feature request, please open an issue on our [GitHub repository](https://github.com/DJ4Earth/DJUICE.jl). If you want to contribute code, please fork the repository and submit a pull request.

## License

DJUICE.jl is released under the MIT License. See the [LICENSE](https://github.com/DJ4Earth/DJUICE.jl/blob/main/LICENSE) file for more details.

## Acknowledgements

This project is inspired by the [ice-sheet and Sea-level System Model (ISSM)](https://issm.jpl.nasa.gov/). Special thanks to the developers of Enzyme.jl for providing an excellent tool for automatic differentiation in Julia.

## Authors:

 - [Mathieu Morlighem](https://github.com/mmorligh), Dartmouth College, USA
 - [Cheng Gong](https://github.com/enigne), Dartmouth College, USA


[build-stable-img]: https://github.com/DJ4Earth/DJUICE.jl/workflows/CI/badge.svg
[build-url]: https://github.com/yourusername/DJUICE.jl/actions








[build-stable-img]: https://github.com/DJ4Earth/DJUICE.jl/workflows/CI/badge.svg
[build-url]: https://github.com/DJ4Earth/DJUICE/actions?query=workflow




