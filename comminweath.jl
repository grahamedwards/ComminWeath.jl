## Simulates U-series comminution system & production of
 # an authigenic "weathered" phase for multiple grain-sizes
 # of multiple ages. ~~ Graham Harper Edwards ~~


## Grain Sizes, Ages, Detrital initials
a = [ 40 10 5 2]' *1e-6 #um -> m
 tcom = [ 10 100 200 300]' * 1000; # ka

AoU = 1. # initial activity ratio of the provenance rock
 AoTh = 1.

## Concentrations.

U_dtr = 300. # ng/g
 U_auth = 1000.# ng/g

## Weathering/Authigenic phase parameters:

k238 = 1e-7 #1/a
 k234 = 1e-7 #weathering rate of each nuclide
 k230 = 5e-8 #negative implies loss, positive addition (auth phase fm)

auth48 = 5 # (234U/238U) of authigenic phase
 auth08 = 1 # (230U/238U) of authigenic phase

## Recoil lengths of pertinent nuclides (nm) in silicates
    #calculated using SRIM 2008 (right?)

L234 = 30 * 1e-9 #nm -> m
 L230 = 37 * 1e-9
 L226 = 34 * 1e-9
## Decay Constants

const l238 = log(2)/ 4.4683e9 # 1/a; Jaffey et al. 1971, Phys. Rev. C
 const l234 = log(2)/ 245620 # 1/a; Cheng et al. 2013, EPSL
 const l230 = log(2)/ 75584 # "  "
 const l232 = log(2)/ 14.05e9 #1/a Audi et al. 1997, Nuc. Phys. A
 const l226 = log(2)/ 1599 # 1/a; Holden, 199, Pure and Applied Chem.

 const SEU = l238/l234
 const SEThU = l238/l230

## Prepare/Pre-allocate for calculaions
nit=10000 # Timesteps for evolution calculation.

A234_ev = zeros(length(a),length(tcom),nit) # grain-size x age result matrix
 A230_ev = copy(A234_ev)

Ndr_238 = U_dtr/238 # mol 238U / g detrital sediment
 Ndr_230i = (AoTh *SEThU) * Ndr_238 # mol230 / g detrital sediment
 Ndr_234i = (AoTh *SEU) * Ndr_238 # mol230 / g detrital sediment

Nau_238 = U_auth/238 # mol 238U / g detrital sediment
 Nau_230 = (auth08 *SEThU) * Nau_238 # mol230 / g detrital sediment
 Nau_234 = (auth48 *SEThU) * Nau_238 # mol230 / g detrital sediment

### --- REMOVE GLOBALS
for z = 1:length(a)

global r = a[z]
global K = 10 # Cartwright, 1962. 10 is ~ best value. K=6 -> sphere
         # 6. + 8. * rand(1)[1] # Randomized grain shapes here for each size class.
global rfr= 7; #surface roughness factor. ~7 for freshly crushed.


# f calculation
 f234 = (L234*K)/(4*2*r)*rfr # from Lee et al. 2010 EPSL
 f230 = (L230*K)/(4*2*r)*rfr

#Numerical Calculation



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



     for j = 2: nit #start on second value.
         N238_ev[j] = N238_ev[j-1] +
                    dt *
                    ( k238 * Nau_238 -
                      l238 * N238_ev[j-1] )
         N234_ev[j] = N234_ev[j-1] +
                    dt *
                    ( k234 * Nau_234 -
                      l234 * N234_ev[j-1] +
                      (1-f234) * l238 * N238_ev[j-1] )
         N230_ev[j] = N230_ev[j-1] +
                    dt *
                    ( k230 * Nau_230 -
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
print(A234)
print(A230)

# Plot Isochrons
#using Pkg; Pkg.add("Plots")
#using Plots; gr() #load Plots and use gr backend.


plot(A230,A234, marker = (:hexagon, 2, 0.6, :black), legend=:topleft)
xlabel!("(230Th/238U)")
ylabel!("(234U/238U)")
# Plot Evolutions
    for i=1:length(a)
        plot!(A230_ev[i,length(tcom),:],A234_ev[i,length(tcom),:])
    end
    current() # to return the plot you produced within the loop
