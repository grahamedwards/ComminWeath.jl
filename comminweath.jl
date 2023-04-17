## Simulates U-series comminution system for multiple grain-sizes
 # of multiple ages, as well as optional
 # production of an authigenic "weathered" phase
 # implantation of U-series nuclides from a high-[U] authigenic rind.
 # ~~ Graham Harper Edwards ~~


include("measured-data.jl")
meas = TaylorIII
include("decay-constants.jl")
include("structs.jl")

cd(@__DIR__)
run_name = "comminwx_out" #comminweath_out
fileout_prefix = "../comminweath_results/"
## Grain Sizes, Ages, Detrital initials:
a = [40. 30. 15. 15.]' *1e-6 #um -> m ~~ Grain diameters
tcom = [ 0 400 1500]' # ka
#Variable [U] with grain size
U_dtr =  [ 210 180 160 130]' #ng/g 
    #|| [U]s correspond with grain sizes in 'a' above
    # should be [300 210 180 160]
    

"""

```julia
function comminweath(diameter, grain, detrital, wxauth, rind; timeseries)
```

Calculates post-comminution U-series evolutions for a grain of a given `diameter` (in μm) over a `timeseries` (`<: AbstractRange`) in years (default=0:1:2e6). 

The history depends on grain morphology as well as the `detrital`, authigenic weathering (`wxauth`), and soluble authigenic `rind` parameters; respectively represented by `mutable struct` types `Grain`, `Detrital`, `WxAuth`, and `Rind`. See the docs of these types for more details.

Outputs a `NamedTuple` containing the `timeseries` and vectors of U elemental and isotopic evolutions with corresponding indices.

| key  | description                                  |
| ---- | -------------------------------------------- |
| `t`  | `timeseries`                                 |
|`A234`| (²³⁴U/²³⁸U) of insoluble (detrital) grain    |
|`A230`| (²³⁰Th/²³⁸U) of insoluble (detrital) grain   |
|`cU`  | [U] of insoluble (detrital) grain            |    
|`etc` | NamedTuple wrapping `rind_` and `bulk_`      |
|      | isotopic and [U] evolutions.                 |
| ---- | -------------------------------------------- |

"""
function comminweath(diameter::Number, grain::Grain, detrital::Detrital, wxauth::WxAuth, rind::Rind; timeseries::AbstractRange=0:1:2e6)
    #d=40 # μm grain size 
    #timeseries = 0:1:2e6

    ## Prepare initial (nano)molar abundances of each isotope in each reservoir.
    Ndr_238 = detrital.cU/238 # mol 238U / g detrital sediment
    Ndr_230i = (detrital.r08 *seThU) * Ndr_238 # mol230 / g detrital sediment, initial
    Ndr_234i = (detrital.r48 *seU) * Ndr_238 # mol230 / g detrital sediment, initial

    Nau_238 = wxauth.cU/238 # mol 238U / g authigenic material
    Nau_230 = (wxauth.r08 *seThU) * Nau_238 # mol230 / g authigenic material
    Nau_234 = (wxauth.r48 *seU) * Nau_238 # mol230 / g authigenic material

    Nr_238 = rind.cU/238 # mol 238U / g nondetrital rind
    Nr_230i = (rind.r08 *seThU) * Nau_238 # mol230 / g nondetrital rind
    Nr_234i = (rind.r48 *seU) * Nau_238 # mol230 / g nondetrital rind


    d, rind_z = (diameter, rind.z) .* 1e-6 # convert μm to m for matching units
    rind_rho = rind.rho * 1e3 # convert kg/m³ to g/m³ for matching units
    L230, L234 = (grain.L230, grain.L234) .* 1e-9 # convert nm to m for matching units.

    # f and S calculation
    f234 = (L234*grain.K)/(4*d)*grain.rfr # from Lee et al. 2010 EPSL
    f230 = (L230*grain.K)/(4*d)*grain.rfr

    S = ( grain.K * grain.rfr ) / ( grain.rho * 1e3 * d ) # m2/g ~~ total surface area
    S_a = ifelse(wxauth.sa_dependent, S, one(S)) # if surface area dependence on,  S_a is an "alteration" S value

    Mr_Md =  rind_rho * S * rind.p * rind_z
    ## Numerical Calculation of U-series evolution
    dt = step(timeseries)

    # pre-allocate detrital evolution vectors
    N238_ev, N234_ev, N230_ev = zero(timeseries), zero(timeseries), zero(timeseries) 
    N238_ev[1], N234_ev[1], N230_ev[1] = Ndr_238, Ndr_234i, Ndr_230i

    # pre-allocate nondetrital rind evolution vectors
    Nr234_ev,Nr230_ev = zero(timeseries), zero(timeseries)
    Nr234_ev[1], Nr230_ev[1] = Nr_234i, Nr_230i
                

    @inbounds for j = 2: length(timeseries) #start on second timestep.

    #### I THINK THE STEP ON `tlin[j]` SHOULD ACTUALLY BE `j-1` ~~ CHANGE AFTER REPRODUCING OUTPUT
        k= wxauth.k * wxauth.k_power * timeseries[j]^(wxauth.k_power-1)
    #### I THINK THE STEP ON `tlin[j]` SHOULD ACTUALLY BE `j-1` ~~ CHANGE AFTER REPRODUCING OUTPUT

        N238_ev[j] = N238_ev[j-1] +
                    dt * (
                    k * S_a * (Nau_238 - N238_ev[j-1]) -    #alteration
                    l238 * N238_ev[j-1] )                 #outgrowth
                                        # see p. 77 in notes on alteration math
        N234_ev[j] = N234_ev[j-1] +
                    dt * (
                    k * S_a * (Nau_234 - N234_ev[j-1]) +    #alteration
                    (L234/4) * l238 * Nr_238 *
                    rind.p * S * rind_rho -               #injection
                    l234 * N234_ev[j-1] +                 #outgrowth
                    (1-f234) * l238 * N238_ev[j-1] )      #ingrowth
                                        # see p. 79 in notes on rind-injection
        N230_ev[j] = N230_ev[j-1] +
                    dt * (
                    k * S_a * (Nau_230 - N230_ev[j-1]) +    #alteration
                    (L230/4) * l234 * Nr234_ev[j-1] *
                    rind.p * S * rind_rho -               #injection
                    l230 * N230_ev[j-1] +                 #outgrowth
                    (1-f230) * l234 * N234_ev[j-1] )      #ingrowth
        # Modified from Eqs. 3 of Cogez et al., 2018

        if rind.evolve
            Nr234_ev[j] = Nr234_ev[j-1] +
                        dt * (
                        -l234 * Nr234_ev[j-1] +                #outgrowth
                        ( 1- 2 * L234 / (4 * rind_z)) * l238 * Nr_238 )
                                                            # ingrowth
            Nr230_ev[j] = Nr230_ev[j-1] +
                        dt * (
                        -l230 * Nr230_ev[j-1] +                #outgrowth
                            (1- 2* L234 / (4* rind_z)) * l238 * Nr234_ev[j-1] )
                                                            # ingrowth
        else
            Nr234_ev[j] = Nr_234i 
            Nr230_ev[j] = Nr_230i
        end
    end

    # Activity Ratios
    A234 = @. (N234_ev / N238_ev) / seU
    A230 = @. (N230_ev / N238_ev) / seThU
    rind_A234 = @. (Nr234_ev / Nr_238) / seU
    rind_A230 = @. (Nr230_ev / Nr_238 ) / seThU
    bulk_A234 = @. ( N234_ev + Nr234_ev * Mr_Md ) / ( N238_ev + Nr_238 * Mr_Md ) / seU
    bulk_A230 = @. ( N230_ev + Nr230_ev * Mr_Md ) / ( N238_ev + Nr_238 * Mr_Md ) / seThU

    # U concentrations
    cU = N238_ev .* 238
    bulk_cU = cU .+ rind.cU * Mr_Md

    (; t=timeseries, A234, A230, cU, etc = (; bulk_A234, bulk_A230, bulk_cU, rind_A234, rind_A230))
end 

"""

```julia
drawdates(dates, evs::NamedTuple)
```

Returns the (²³⁴U/²³⁸U) (`A234`), (²³⁰Th/²³⁸U) (`A230`), and U concentration (`cU`) of a sediment grain for a collection of `dates` in a U-series history (`evs`) from the `comminweath` function.
"""
function drawdates(dates,evs::NamedTuple)
    n = length(dates)
    d = dates .* 1e3

    A234, A230, cU = Vector{eltype(evs.A234)}(undef,n), Vector{eltype(evs.A230)}(undef,n), Vector{eltype(evs.cU)}(undef,n)

    @inbounds for i in eachindex(d)
        @assert d[i] ∈ evs.t "$(dates[i]) ka is not a simulated timestep. Re-run comminweath with an expanded timeseries."
        t = searchsortedfirst(evs.t,d[i])
        A234[i], A230[i], cU[i] = evs.A234[t], evs.A230[t], evs.cU[t]
    end

    (; A234, A230, cU)
end

include("regressions.jl")
include("visualize.jl")
