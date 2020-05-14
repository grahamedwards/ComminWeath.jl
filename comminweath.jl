## Simulates U-series comminution system & production of
 # an authigenic "weathered" phase for multiple grain-sizes
 # of multiple ages. ~~ Graham Harper Edwards ~~

#import DataFrames; import GLM; using Plots; gr()
run_name = "comminweath_out"
fileout_prefix = "~/Dropbox/GHE/Science/4_U-Series/Comminution/Comminution_Glacial/Taylor Valley/model_results/"
## Grain Sizes, Ages, Detrital initials:
a = [ 40 25 15 ]' *1e-6 #um -> m ~~ Grain diameters
 tcom = [10 100 300 500 1000]' * 1000; # ka

AoU = 0.95 # initial activity ratio of the provenance rock
 AoTh = 0.95

## Concentrations:

U_dtr = 200. # ng/g
 U_auth = 1000.# ng/g

## Weathering/Authigenic phase parameters:

k_coeff = 5e-6 #1e-6 # g/[m2 a] ~~ Rate of alteration to authigenic phase
 k_power = 1/2
 auth48 = 10 # (234U/238U) of authigenic phase
 auth08 = 10 # (230U/238U) of authigenic phase

## Grain and alpha-recoil physical parameters
K = 10 # Cartwright, 1962. 10 is ~ best value. K=6 -> sphere
         # 6. + 8. * rand(1)[1] # Randomized grain shapes here for each size class.
rfr = 7 #surface roughness factor. ~7 for freshly crushed.
rho = 2650 * 1e3 #kg/m3 --> g/m3 ~~ sediment density

#Recoil lengths of pertinent nuclides (nm) in silicates
    #calculated using SRIM 2008 (right?)
L234 = 30 * 1e-9 #nm -> m
 L230 = 37 * 1e-9
 L226 = 34 * 1e-9
## Decay Constants

const l238 = log(2)/ 4.4683e9 # 1/a; Jaffey et al. 1971, Phys. Rev. C
 const l234 = log(2)/ 245620 # 1/a; Cheng et al. 2013, EPSL
 const l230 = log(2)/ 75584 # "  "
 const l232 = log(2)/ 14.05e9 #1/a Audi et al. 1997, Nuc. Phys. A
 const l226 = log(2)/ 1599 # 1/a; Holden, 1990, Pure and Applied Chem.

 const SEU = l238/l234
 const SEThU = l238/l230

## Prepare/Pre-allocate for calculaions
nit=10000 # Timesteps for evolution calculation.

A234_ev = zeros(length(a),length(tcom),nit) # grain-size x age result matrix
 A230_ev = copy(A234_ev)

Ndr_238 = U_dtr/238 # mol 238U / g detrital sediment
 Ndr_230i = (AoTh *SEThU) * Ndr_238 # mol230 / g detrital sediment, initial
 Ndr_234i = (AoU *SEU) * Ndr_238 # mol230 / g detrital sediment, initial

Nau_238 = U_auth/238 # mol 238U / g authigenic material
 Nau_230 = (auth08 *SEThU) * Nau_238 # mol230 / g authigenic material
 Nau_234 = (auth48 *SEU) * Nau_238 # mol230 / g authigenic material

## Calculate U-series evolution for each grain size
for z = 1:length(a)

d = a[z]
# f and S calculation
 f234 = (L234*K)/(4*d)*rfr # from Lee et al. 2010 EPSL
 f230 = (L230*K)/(4*d)*rfr

 S = ( K * rfr ) / ( rho * d ) # m2/g ~~ total surface area

## Numerical Calculation of U-series evolution

for i = 1:length(tcom)
     t=tcom[i]
     dt = t/(nit-1)
     tlin=LinRange(0,t,nit)
     N238_ev = zeros(nit,1)
     N234_ev = copy(N238_ev)
     N230_ev = copy(N238_ev)
     N238_ev[1] = Ndr_238
     N234_ev[1] = Ndr_234i
     N230_ev[1] = Ndr_230i



     for j = 2: nit #start on second timestep.

         k= k_coeff * k_power * tlin[j]^(k_power-1)

         N238_ev[j] = N238_ev[j-1] +
                    dt *
                    ( k * S * (Nau_238 - N238_ev[j-1]) -  #see p. 77 in notes
                      l238 * N238_ev[j-1] )

         N234_ev[j] = N234_ev[j-1] +
                    dt *
                    ( k * S * (Nau_234 - N234_ev[j-1]) -
                      l234 * N234_ev[j-1] +
                      (1-f234) * l238 * N238_ev[j-1] )

         N230_ev[j] = N230_ev[j-1] +
                    dt *
                    ( k * S * (Nau_230 - N230_ev[j-1]) -
                      l230 * N230_ev[j-1] +
                      (1-f230) * l234 * N234_ev[j-1] )
         # Modified from Eqs. 3 of Cogez et al., 2018
     end
     A234_ev[z,i,:] = (N234_ev ./ N238_ev) ./ SEU
     A230_ev[z,i,:] = (N230_ev ./ N238_ev) ./ SEThU
     #z~grain size, i~age
end
end
A230 = copy(A230_ev[:,:,end])
A234 = copy(A234_ev[:,:,end])


## Calculate Slopes and r2 of each

r2s = zeros(length(tcom),1)
slope_int = zeros(length(tcom),2)
for i = 1: length(A234[1,:])
    df = DataFrames.DataFrame(x=A230[:,i],y=A234[:,i])
    reg = GLM.lm(@formula(y~x),df)
    r2s[i] = r2(reg)
    slope_int[i,:] = coef(reg) # col1: intercept, col2: slope
end
# Plot Isochrons
#using Plots; gr() #load Plots with gr backend.
p1 = plot(A230,A234, marker = (:hexagon, 2, 0.6, :black),
                     legend=:bottomright, label = Int.(tcom'/1000),
                     legendtitle="Age (ka)")
xaxis!("(230Th/238U)")
ylabel!("(234U/238U)")
# Plot Evolutions
    for i=1:length(a)
        p1=plot!(A230_ev[i,length(tcom),:],A234_ev[i,length(tcom),:],
        color= :gray, label="")
    end
p2 = plot(tcom/1000,slope_int[:,2], label ="",ylim=(-1,1))

ylabel!("slope");xlabel!("Age (ka)")

plot(p1,p2, layout=(1,2))
current()
savefig(string(fileout_prefix,run_name,".pdf"))

print(r2s)
