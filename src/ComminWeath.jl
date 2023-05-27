module ComminWeath

include("decay-constants.jl")

export Grain, Rind, WxAuth, Detrital
include("structs.jl")

export comminweath, drawdate, drawdates, linreg, calcslopes, calcU
include("maths.jl")

export Tanaka2015, TaylorI, TaylorIII, TaylorIV
include("measured-data.jl")

# VisualizeComminWeath is a separate module, which can probably be exposed with:  
push!(LOAD_PATH, @__DIR__)

end
