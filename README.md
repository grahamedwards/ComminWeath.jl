# ComminWeath.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GrahamEdwards.github.io/ComminWeath.jl/dev/)
[![Build Status](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl)

Exploring the effects of comminution and weathering on the U-series systematics of sediments.

This package contains a suite of functions that simulate U-series comminution systematics of fine-grained sediments. Incorporates replacement of detrital material with an insoluble authigenic "weathered" phase and implantation of U-series nuclides from a soluble authigenic rind.

## Usage

The easiest way to use/explore the `ComminWeath.jl` package will be through the online binder jupyter notebook, which is currently under construction. 

Alternatively, the online [Documentation](https://grahamedwards.github.io/ComminWeath.jl/dev/) details all of the contents of the package and their functionality.

## Installation

To install, just type into the Julia REPL `]` to enter the built-in package manager and then type `add https://github.com/grahamedwards/ComminWeath.jl` and hit enter.

After installing, just type `using ComminWeath` to use. 

### Visualizations
The visualization functions in the `examples/` directory package rely on [Makie.jl](https://docs.makie.org/stable/) to make visualizations. If you wish to use these, you'll need to install a Makie backend in your Julia environment. Makie has a couple of different backends available (e.g. [CairoMakie](https://docs.makie.org/stable/documentation/backends/cairomakie/), [GLMakie](https://docs.makie.org/stable/documentation/backends/glmakie/)). I recommend `CairoMakie` for our purposes. To install, just type into the Julia REPL `]add CairoMakie` and hit enter. It may take a little bit of time to install.

Then, before using any of the , type `using CairoMakie` into the REPL and hit enter.

---

Please reach out to me if you have any questions or encounter any bugs!