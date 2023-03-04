## Simulates U-series comminution system for multiple grain-sizes
 # of multiple ages, as well as optional
 # production of an authigenic "weathered" phase
 # implantation of U-series nuclides from a high-[U] authigenic rind.
 # ~~ Graham Harper Edwards ~~


include("measured-data.jl")
meas = TaylorIII

run_name = "comminwx_out" #comminweath_out
fileout_prefix = "~/Dropbox/GHE/Science/4_U-Series/Comminution/Comminution_Glacial/Taylor Valley/model_results/"
## Grain Sizes, Ages, Detrital initials:
a = [40. 30. 15. 15.]' *1e-6 #um -> m ~~ Grain diameters
 tcom = [ 0 400 1500]' * 1000; # ka

AoU = 1. # initial (234U/238U) of the provenance rock
 AoTh = 1. # initial (230Th/238U) of the provenance rock

## Concentrations:

#Assume fixed [U] in detrital component
#  U_dtr = 300. * ones(length(a))# ng/g || fixed detrital [U]
#Variable [U] with grain size
  U_dtr =  [ 210 180 160 130]' #ng/g
        #|| [U]s correspond with grain sizes in 'a' above
        # should be [300 210 180 160]
# Assume fixed [U] in authigenic component
 U_auth = 1000. # ng/g #1000 usually
 U_rind = 2000. # ng/g

## Weathering/Authigenic phase parameters:
wxauth = (
    sa_dependent = true, #Surface area dependent? ON: 1 || OFF: 0
    k = 0.0e-6, # ~~ Rate of alteration to authigenic phase (g/[m2 a] OR g/[g a]), (initial rate if a power-law relationship is used)
    k_power = 1, # Exponent on k, if using a power-law weathering rate. Not explored in this study, so set to 1 for constant weathering rate.
    r48 = 2, # (²³⁴U/²³⁸U) of authigenic phase
    r08 = 3 # (²³⁰Th/²³⁸U) of authigenic phase
    )

## Authigenic rind parameters
rind = (
    evolve = true, # Does the rind evolve?
    p = .5,  # proportion of grain covered by rind (q on notes p.79)
    z = .1 * 1e-6, # um -> m ~~ thickness of nondetrital rind.
    rho = 2500 * 1e3, #kg/m3 -> g/m3 ~~ nondetrital rind phase density
    r48 = 3  # (234U/238U) of nondetrital rind
    r08 = 4.6# (230Th/238U) of nondetrital rind
)
## Grain and alpha-recoil physical parameters
grain = (
    K = 10, # Cartwright, 1962. 10 is ~ best value. K=6 -> sphere
    rfr = 7,  #surface roughness factor. ~7 for freshly crushed.
    rho = 2650 * 1e3, #kg/m3 -> g/m3 ~~ sediment density
    L234 = 34 * 1e-9, # nm -> m ~~ recoil length of ²³⁴U in framework silicates (SRIM in 2008 -> 30 nm)
    L230 = 37 * 1e-9 # nm -> m ~~ recoil length of ²³⁰Th in framework silicates (SRIM in 2008 -> 37 nm)
    )

## Prepare/Pre-allocate for calculaions
nit=100000 # Timesteps for evolution calculation.

A234_ev = zeros(length(a),length(tcom),nit) # grain-size x age result matrix
 A230_ev = copy(A234_ev)
 bulk_A234_ev = copy(A234_ev)
 bulk_A230_ev = copy(A234_ev)
 cUd_ev = zeros(length(a),length(tcom),nit)
 cU_bulk = zeros(length(a), length(tcom))
 r_A234_ev = zeros(length(a),length(tcom),nit)
 r_A230_ev = zeros(length(a),length(tcom),nit)
Ndr_238 = U_dtr/238 # mol 238U / g detrital sediment
 Ndr_230i = (AoTh *SEThU) * Ndr_238 # mol230 / g detrital sediment, initial
 Ndr_234i = (AoU *SEU) * Ndr_238 # mol230 / g detrital sediment, initial

Nau_238 = U_auth/238 # mol 238U / g authigenic material
 Nau_230 = (wxauth.r08 *SEThU) * Nau_238 # mol230 / g authigenic material
 Nau_234 = (wxauth.r48 *SEU) * Nau_238 # mol230 / g authigenic material

 Nr_238 = U_rind/238 # mol 238U / g nondetrital rind
  Nr_230i = (rind.r08 *SEThU) * Nau_238 # mol230 / g nondetrital rind
  Nr_234i = (rind.r48 *SEU) * Nau_238 # mol230 / g nondetrital rind

## Calculate U-series evolution for each grain size
for z = 1:length(a)

d = a[z]
# f and S calculation
 f234 = (L234*grain.K)/(4*d)*rfr # from Lee et al. 2010 EPSL
 f230 = (L230*grain.K)/(4*d)*rfr

 S = ( grain.K * grain.rfr ) / ( grain.rho * d ) # m2/g ~~ total surface area
        if wxauth.sa_dependent  #surface area dependence on
            S_a = S  #S_a is an "alteration" S value
        else #surface area dependence off
            S_a = 1
        end
 Mr_Md =  rind.rho * S * rind.p * rind.z
## Numerical Calculation of U-series evolution

    for i = 1:length(tcom)
        t=tcom[i]
        dt = t/(nit-1)
        tlin=LinRange(0,t,nit)
        N238_ev = zeros(nit,1) ; N238_ev[1] = Ndr_238[z]
        N234_ev = zeros(nit,1) ; N234_ev[1] = Ndr_234i[z]
        N230_ev = zeros(nit,1) ; N230_ev[1] = Ndr_230i[z]
            #pre-allocate detrital evolution vectors
        #Nr238_ev = zeros(nit,1) ; Nr238_ev[1] = Nr_238 # assuming does not change
        Nr234_ev = zeros(nit,1) ; Nr234_ev[1] = Nr_234i
        Nr230_ev = zeros(nit,1) ; Nr230_ev[1] = Nr_230i
            #pre-allocate nondetrital rind evolution vectors

        for j = 2: nit #start on second timestep.

            k= wxauth.k * k_power * tlin[j]^(k_power-1)

            N238_ev[j] = N238_ev[j-1] +
                        dt * (
                        k * S_a * (Nau_238 - N238_ev[j-1]) -    #alteration
                        l238 * N238_ev[j-1] )                 #outgrowth
                                            # see p. 77 in notes on alteration math
            N234_ev[j] = N234_ev[j-1] +
                        dt * (
                        k * S_a * (Nau_234 - N234_ev[j-1]) +    #alteration
                        (L234/4) * l238 * Nr_238 *
                        rind.p * S * rind.rho -               #injection
                        l234 * N234_ev[j-1] +                 #outgrowth
                        (1-f234) * l238 * N238_ev[j-1] )      #ingrowth
                                            # see p. 79 in notes on rind-injection
            N230_ev[j] = N230_ev[j-1] +
                        dt * (
                        k * S_a * (Nau_230 - N230_ev[j-1]) +    #alteration
                        (L230/4) * l234 * Nr234_ev[j-1] *
                        rind.p * S * rind.rho -               #injection
                        l230 * N230_ev[j-1] +                 #outgrowth
                        (1-f230) * l234 * N234_ev[j-1] )      #ingrowth
            # Modified from Eqs. 3 of Cogez et al., 2018

            if rind.evolve
                Nr234_ev[j] = Nr234_ev[j-1] +
                            dt * (
                            -l234 * Nr234_ev[j-1] +                #outgrowth
                            ( 1- 2 * L234 / (4 * rind.z)) * l238 * Nr_238 )
                                                                # ingrowth
                Nr230_ev[j] = Nr230_ev[j-1] +
                            dt * (
                            -l230 * Nr230_ev[j-1] +                #outgrowth
                                (1- 2* L234 / (4* rind.z)) * l238 * Nr234_ev[j-1] )
                                                                # ingrowth
            else
                Nr234_ev[j] = Nr_234i ; Nr230_ev[j] = Nr_230i
            end
        end

        A234_ev[z,i,:] = (N234_ev ./ N238_ev) ./ SEU
        A230_ev[z,i,:] = (N230_ev ./ N238_ev) ./ SEThU
        cUd_ev[z,i,:] = N238_ev*238

        r_A234_ev[z,i,:] = (Nr234_ev ./ Nr_238) ./ SEU
        r_A230_ev[z,i,:] = (Nr230_ev ./ Nr_238 ) ./ SEThU

        bulk_A234_ev[z,i,:] = @. ( N234_ev + Nr234_ev * Mr_Md ) /
                                    ( N238_ev + Nr_238 * Mr_Md ) / SEU
        bulk_A230_ev[z,i,:] = @. ( N230_ev + Nr230_ev * Mr_Md ) /
                                    ( N238_ev + Nr_238 * Mr_Md ) / SEThU
        #z~grain size, i~age
    end
    cU_bulk[z,:] = cUd_ev[z,:,end] .+ ( U_rind * Mr_Md )
end

A230 = copy(A230_ev[:,:,end])
A234 = copy(A234_ev[:,:,end])
cUd = copy(cUd_ev[:,:,end])

r_A234 = copy(r_A234_ev[:,:,end])
r_A230 = copy(r_A230_ev[:,:,end])

bulk_A230 = copy(bulk_A230_ev[:,:,end])
bulk_A234 = copy(bulk_A234_ev[:,:,end])


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

## Plot Isochrons & Age Calculations
using Plots; gr() #load Plots with gr backend.
p1 = plot(A230,A234,
                     legend=false, label = Int.(tcom'/1000),
                     legendtitle="Age (ka)", marker = (:hexagon, 2, 0.6, :black))
xaxis!("(230Th/238U)")
ylabel!("(234U/238U)")
# Plot Evolutions
    for i=1:length(a)
        global p1 =plot!(A230_ev[i,length(tcom),:],A234_ev[i,length(tcom),:],
        color= :gray, linewidth = 0.1, label="")
    end
p1=plot!(meas.r08,meas.r48,xerr=meas.u08, yerr=meas.u48)
# Plot slope vs. age
p2 = plot(tcom/1000,slope_int[:,2], label ="") #, ylim=(0,0.6))
ylabel!("slope");xlabel!("Age (ka)")

# Plot [U]_detrital vs. grain size
p3 = plot(repeat(a,1,length(tcom))*1e6,cUd/1e3,
legend=false, label = Int.(tcom'/1000),
legendtitle="Age (ka)", marker = (:hexagon, 2, 0.6, :black))
ylabel!("[U] (ppm)");xlabel!("Grain size (um)")

# Calculate fraction clay from cUd
auth_pct = 100*(cUd .- repeat(U_dtr,1,length(tcom))) ./ (U_auth .- U_dtr)
p4 = plot(repeat(a,1,length(tcom))*1e6,auth_pct,
legend=false, label = Int.(tcom'/1000),
legendtitle="Age (ka)",marker = (:hexagon, 2, 0.6, :black))
ylabel!("Percent authigenic phase (%)");xlabel!("Grain size (um)")

plot_out=plot(p1,p2,p3,p4, layout=(2,2))
current()
savefig(string(fileout_prefix,run_name,".pdf"))
display(plot_out)
