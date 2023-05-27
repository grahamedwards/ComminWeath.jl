module ComminWeath

include("decay-constants.jl")

export Grain, Rind, WxAuth, Detrital
include("structs.jl")

export comminweath, drawdate, drawdates, linreg, calcslopes, calcU
include("maths.jl")

export Tanaka2015, TaylorI, TaylorIII, TaylorIV
include("measured-data.jl")

import ColorSchemes # for visualizations
using Makie

export plotUseries, plotslopes, plotcU, plotauthreplace
include("visualize.jl")


end