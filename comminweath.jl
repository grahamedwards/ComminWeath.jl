## Simulates U-series comminution system for multiple grain-sizes
 # of multiple ages, as well as optional
 # production of an authigenic "weathered" phase
 # implantation of U-series nuclides from a high-[U] authigenic rind.
 # ~~ Graham Harper Edwards ~~


include("measured-data.jl")
meas = TaylorIII
include("decay-constants.jl")

cd(@__DIR__)
run_name = "comminwx_out" #comminweath_out
fileout_prefix = "../comminweath_results/"
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
    sa_dependent = true, # Surface area dependent weathering
    k = 0.0e-6, # Rate of alteration to authigenic phase (g/[m2 a] OR g/[g a]), (initial rate if a power-law relationship is used)
    k_power = 1, # Exponent on k, if using a power-law weathering rate. Not explored in this study, so set to 1 for constant weathering rate, but there if you want it!
    r48 = 2, # (²³⁴U/²³⁸U) of authigenic phase
    r08 = 3 # (²³⁰Th/²³⁸U) of authigenic phase
    )

## Authigenic rind parameters
rind = (
    evolve = true, # Does the rind evolve?
    p = .5,  # proportion of grain covered by rind (q on notes p.79)
    z = .1 * 1e-6, # um -> m ~~ thickness of nondetrital rind.
    rho = 2500 * 1e3, #kg/m3 -> g/m3 ~~ nondetrital rind phase density
    r48 = 3,  # (234U/238U) of nondetrital rind
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
for z = eachindex(a)

d = a[z]
# f and S calculation
 f234 = (grain.L234*grain.K)/(4*d)*grain.rfr # from Lee et al. 2010 EPSL
 f230 = (grain.L230*grain.K)/(4*d)*grain.rfr

 S = ( grain.K * grain.rfr ) / ( grain.rho * d ) # m2/g ~~ total surface area
        if wxauth.sa_dependent  #surface area dependence on
            S_a = S  #S_a is an "alteration" S value
        else #surface area dependence off
            S_a = 1
        end
 Mr_Md =  rind.rho * S * rind.p * rind.z
## Numerical Calculation of U-series evolution

    for i = eachindex(tcom)
        t=tcom[i]
        tlin=LinRange(0,t,nit)
        dt = step(tlin)
        N238_ev = zeros(nit,1) ; N238_ev[1] = Ndr_238[z]
        N234_ev = zeros(nit,1) ; N234_ev[1] = Ndr_234i[z]
        N230_ev = zeros(nit,1) ; N230_ev[1] = Ndr_230i[z]
            #pre-allocate detrital evolution vectors
        #Nr238_ev = zeros(nit,1) ; Nr238_ev[1] = Nr_238 # assuming does not change
        Nr234_ev = zeros(nit,1) ; Nr234_ev[1] = Nr_234i
        Nr230_ev = zeros(nit,1) ; Nr230_ev[1] = Nr_230i
            #pre-allocate nondetrital rind evolution vectors

        for j = 2: nit #start on second timestep.

            k= wxauth.k * wxauth.k_power * tlin[j]^(wxauth.k_power-1)

            N238_ev[j] = N238_ev[j-1] +
                        dt * (
                        k * S_a * (Nau_238 - N238_ev[j-1]) -    #alteration
                        l238 * N238_ev[j-1] )                 #outgrowth
                                            # see p. 77 in notes on alteration math
            N234_ev[j] = N234_ev[j-1] +
                        dt * (
                        k * S_a * (Nau_234 - N234_ev[j-1]) +    #alteration
                        (grain.L234/4) * l238 * Nr_238 *
                        rind.p * S * rind.rho -               #injection
                        l234 * N234_ev[j-1] +                 #outgrowth
                        (1-f234) * l238 * N238_ev[j-1] )      #ingrowth
                                            # see p. 79 in notes on rind-injection
            N230_ev[j] = N230_ev[j-1] +
                        dt * (
                        k * S_a * (Nau_230 - N230_ev[j-1]) +    #alteration
                        (grain.L230/4) * l234 * Nr234_ev[j-1] *
                        rind.p * S * rind.rho -               #injection
                        l230 * N230_ev[j-1] +                 #outgrowth
                        (1-f230) * l234 * N234_ev[j-1] )      #ingrowth
            # Modified from Eqs. 3 of Cogez et al., 2018

            if rind.evolve
                Nr234_ev[j] = Nr234_ev[j-1] +
                            dt * (
                            -l234 * Nr234_ev[j-1] +                #outgrowth
                            ( 1- 2 * grain.L234 / (4 * rind.z)) * l238 * Nr_238 )
                                                                # ingrowth
                Nr230_ev[j] = Nr230_ev[j-1] +
                            dt * (
                            -l230 * Nr230_ev[j-1] +                #outgrowth
                                (1- 2* grain.L234 / (4* rind.z)) * l238 * Nr234_ev[j-1] )
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

include("regressions.jl")
include("visualize.jl")
