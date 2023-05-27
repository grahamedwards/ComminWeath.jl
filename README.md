# ComminWeath.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GrahamEdwards.github.io/ComminWeath.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GrahamEdwards.github.io/ComminWeath.jl/dev/)
[![Build Status](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl)

Exploring the effects of comminution and weathering on the U-series systematics of sediments.

This package contains a suite of functions that simulate U-series comminution systematics of fine-grained sediments. Incorporates replacement of detrital material with an insoluble authigenic "weathered" phase and implantation of U-series nuclides from a soluble authigenic rind.  


## Visualizations
This package relies on [Makie.jl](https://docs.makie.org/stable/) to make its visualizations. Since plotting packages take a bit of computational overhead and slow things down, I'm leaving it up to you to install a Makie backend in your Julia environment. Makie has a couple of different backends available (e.g. [CairoMakie](https://docs.makie.org/stable/documentation/backends/cairomakie/), [GLMakie](https://docs.makie.org/stable/documentation/backends/glmakie/)). I recommend `CairoMakie` for our purposes. To install, just type into the Julia REPL `]add CairoMakie` and hit enter. It may take a little bit of time to install.

Then, before making your plot, type `using CairoMakie` into the REPL and hit enter.