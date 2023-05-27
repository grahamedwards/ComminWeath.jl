# ComminWeath.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GrahamEdwards.github.io/ComminWeath.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GrahamEdwards.github.io/ComminWeath.jl/dev/)
[![Build Status](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl)

Exploring the effects of comminution and weathering on the U-series systematics of sediments.

This package contains a suite of functions that simulate U-series comminution systematics of fine-grained sediments. Incorporates replacement of detrital material with an insoluble authigenic "weathered" phase and implantation of U-series nuclides from a soluble authigenic rind.  


## Visualizations
To use any of the visualizations, you will need to install a version of Makie.jl in your Julia environment. Makie has a couple of different backends available. I recommend `CairoMakie` for our purposes. To install, just type into the Julia REPL `]add CairoMakie` and hit enter. It may take a while to install.

Then, before making your plot, type `using CairoMakie` into the REPL and hit enter.