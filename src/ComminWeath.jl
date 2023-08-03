module ComminWeath

include("decay-constants.jl")

export Grain, Rind, WxAuth, Detrital
include("structs.jl")

export comminweath, drawdate, drawdates, linreg, calcslopes, calcU
include("maths.jl")

export TaylorI, TaylorIII, TaylorIV
include("measured-data.jl")

end