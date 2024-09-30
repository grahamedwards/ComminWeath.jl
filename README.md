# ComminWeath.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GrahamEdwards.github.io/ComminWeath.jl/dev/)
[![Build Status](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GrahamEdwards/ComminWeath.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/GrahamEdwards/ComminWeath.jl)

Exploring the effects of comminution and weathering on the U-series systematics of sediments.

This package contains functions that simulate U-series comminution systematics of fine-grained sediments. The models incorporate replacement of detrital material with an insoluble authigenic "weathered" phase and implantation of U-series nuclides from a soluble authigenic rind.

This package accompanies a forthcoming manuscript... stay tuned!

## Usage

The easiest way to use/explore the `ComminWeath.jl` package is through the `examples.ipynb` notebook (in the [/examples](https://github.com/grahamedwards/ComminWeath.jl/tree/main/examples) directory), which runs on a JupyterHub server hosted by the [Binder project](https://mybinder.org/): 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/grahamedwards/ComminWeath.jl/main?labpath=examples%2Fexamples.ipynb) &larr;click this button to load

⚠️ *This may be slow to load or require a few tries, please be patient* ⚠️

Alternatively, the online [Documentation](https://grahamedwards.github.io/ComminWeath.jl/dev/) details all of the contents of the package and their functionality.

## Installation

To install `ComminWeath.jl` on your own computer, just type `]` into the Julia REPL to enter the built-in package manager and then type
```
add https://github.com/grahamedwards/ComminWeath.jl
```
and hit enter.


After installing, just type `using ComminWeath` to use. 

## Visualizations

The easiest way to explore these codes is through the previously mentioned [Binder notebook](https://mybinder.org/v2/gh/grahamedwards/ComminWeath.jl/main?labpath=examples%2Fexamples.ipynb).

Alternatively, the [examples/](examples/) folder contains a Jupyter notebook of examples (`examples.ipynb`) and accompanying suite of visualization functions (`visualizations.jl`). Downloading this directory or cloning the repository will allow you to run that notebook on your own computer. 

### A note for plotting on your own computer:
The visualization functions rely on the [Makie.jl](https://docs.makie.org/stable/) plotting package. If you wish to use any of these included plotting functions, you will need to install a Makie backend in your Julia environment. Makie has a couple of different backends available (e.g. [CairoMakie](https://docs.makie.org/stable/documentation/backends/cairomakie/), [GLMakie](https://docs.makie.org/stable/documentation/backends/glmakie/)). 

I recommend `CairoMakie` for our purposes and include that as a package dependency for use in [the Binder notebook](https://mybinder.org/v2/gh/grahamedwards/ComminWeath.jl/main?labpath=examples%2Fexamples.ipynb). To install, just type `]` into the Julia REPL to enter the built-in package manager and then type `add CairoMakie` (or your backend of choice) and hit enter. It may take a little bit of time to install. Then, before using any of the plotting functions, type `using CairoMakie` into the REPL (or notebook cell) and hit enter.

---

Please reach out to me if you have any questions or encounter any bugs!
