## Calculate Slopes and r2 of each
import DataFrames, GLM
using StatsModels
r2s = zeros(length(tcom),1)
slope_int = zeros(length(tcom),2)
for i = 1: length(A234[1,:])
    df = DataFrames.DataFrame(x=A230[:,i],y=A234[:,i])
    reg = GLM.lm(@formula(y~x),df)
    r2s[i] = GLM.r2(reg)
    slope_int[i,:] = GLM.coef(reg) # col1: intercept, col2: slope
end

print(r2s)