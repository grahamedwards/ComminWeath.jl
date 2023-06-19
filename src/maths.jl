# comminweath, drawdates    

"""

```julia
function comminweath(diameter, grain, detrital, wxauth, rind; timeseries)
```

Calculates post-comminution U-series evolutions for a grain of a given `diameter` (in μm) over a `timeseries` (`<: AbstractRange`) in years (default=0:1:2e6). 

The history depends on grain morphology as well as the `detrital`, authigenic weathering (`wxauth`), and soluble authigenic `rind` parameters; respectively represented by `mutable struct` types `Grain`, `Detrital`, `WxAuth`, and `Rind`. See the docs of these types for more details.

Outputs a `NamedTuple` containing the `timeseries` and vectors of U elemental and isotopic evolutions with corresponding indices.

| key  | description                                                          |
| :--- | :------------------------------------------------------------------- |
| `t`  | `timeseries`                                                         |
|`A234`| (²³⁴U/²³⁸U) of insoluble (detrital) grain                            |
|`A230`| (²³⁰Th/²³⁸U) of insoluble (detrital) grain                           |
|`cU`  | [U] of insoluble (detrital) grain                                    |    
|`etc` | NamedTuple wrapping `rind_` and `bulk_` isotopic and [U] evolutions. |


"""
function comminweath(diameter::Number, grain::Grain, detrital::Detrital, wxauth::WxAuth, rind::Rind; timeseries::AbstractRange=0:1:2e6)
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
    dt = float(step(timeseries))

    # pre-allocate detrital and non-detrital rind evolution vectors
    N238_ev = zeros(float(eltype(timeseries)),length(timeseries))
    N234_ev = copy(N238_ev)
    N230_ev = copy(N238_ev)
    Nr234_ev = copy(N238_ev)
    Nr230_ev = copy(N238_ev)

    N238_ev[1], N234_ev[1], N230_ev[1], Nr234_ev[1], Nr230_ev[1] = Ndr_238, Ndr_234i, Ndr_230i, Nr_234i, Nr_230i

    @inbounds for j = 2: length(timeseries) #start on second timestep.

        k= wxauth.k * wxauth.k_power * timeseries[j-1]^(wxauth.k_power-1)

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

Returns the (²³⁴U/²³⁸U) (`A234`), (²³⁰Th/²³⁸U) (`A230`), and U concentration (`cU`) of a sediment grain for a collection of `dates` (in ka), given a U-series history (`evs`) from the `comminweath` function.
"""
function drawdates(dates::AbstractVector,evs::NamedTuple)
    n = length(dates)
    
    A234, A230, cU = similar(evs.A234,n), similar(evs.A230,n), similar(evs.cU,n)

    @inbounds for i in eachindex(dates)
        A234[i], A230[i], cU[i] = drawdate(dates[i],evs)
    end
    (; A234, A230, cU)
end

"""

```julia
drawdates(dates, evs::NamedTuple)
```

Returns the (²³⁴U/²³⁸U) (`A234`), (²³⁰Th/²³⁸U) (`A230`), and U concentration (`cU`) of a sediment grain for a `dates` (in ka), given a U-series history (`evs`) from the `comminweath` function.

"""
function drawdate(d::Number,evs::NamedTuple)
    d *= 1e3
    @assert d ∈ evs.t "$(d*1e-3) ka is not a simulated timestep. Re-run comminweath with an expanded timeseries."
    
    t = searchsortedfirst(evs.t,d)
    (; A234=evs.A234[t], A230=evs.A230[t], cU=evs.cU[t])
end

"""

```julia
linreg(x::Vector,y::Vector)
```

Calculate a linear regression by QR decomposition for values of independent variable `x` and dependent variable `y`.

Returns a `NamedTuple` with slope `m`, intercept `b`, and coefficient of determination `r²`
"""
function linreg(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y)
    # Calculate linear regression by QR decomposition
    b,m = [ones(length(x)) x] \ y # as in `y = mx+b`

    # Calculate r²
    sst = ssr = zero(m)
    μy = sum(y)/length(y) 
    
    @inbounds for  i in eachindex(y)
        ssr += (y[i] - (m*x[i] + b))^2
        sst += (y[i] - μy)^2
    end
    r² = 1 - ssr/sst

    (; m,b,r²)
end 


"""

```julia
calcslopes(t; a, cU, g=Grain(), d=Detrital(), wx=WxAuth(), r=Rind())
```

Calculate the (²³⁴U/²³⁸U)-(²³⁰Th/²³⁸U) slopes of grain sizes in `a` (default= `[10.,20.,30.,40.]`) at times given in `t` (ka).

A Vector of U concentrations (ng/g) may be passed to `cU` to simulate grain-size dependent U concentrations, but this requires `length(cU) == length(a)`.

Custom parameterizations for Grain, Detrital, WxAuth, and Rind may be passed to `g`, `d`, `wx` and `r`, respectively.

see also: [`Grain`](@ref), [`Detrital`](@ref), [`WxAuth`](@ref), [`Rind`](@ref)

"""
function calcslopes(t::AbstractVector; a::Vector=[10.,20.,30.,40.],cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    isempty(cU) || @assert length(cU) == length(a)
    A234 = zeros(Float64, length(a),length(t))
    A230 = copy(A234)
    cwouts = Vector{NamedTuple}(undef,length(a))
    slopes = similar(t,Float64)

    for i in eachindex(a)
        isempty(cU) || (d.cU=cU[i])
        cwouts[i] = comminweath(a[i],g,d,wx,r)
        for j in eachindex(t)
            Uout = drawdate(t[j],cwouts[i])
            A234[i,j] = Uout.A234
            A230[i,j] = Uout.A230
        end
    end
    for i in eachindex(t)
        x = view(A230,:,i)
        y = view(A234,:,i)
        if !iszero(t[i])
            slopes[i] = linreg(x,y).m
        end
    end
    slopes
end

"""

```julia
calcU(a, t; cU, g=Grain(), d=Detrital(), wx=WxAuth(),r=Rind())
```

Calculate the U abundance for sediments with grain sizes in `a` at times in `t` (ka).

A Vector of U concentrations (ng/g) may be passed to `cU` to simulate grain-size dependent U concentrations, but this requires `length(cU) == length(a)`.

Custom parameterizations for Grain, Detrital, WxAuth, and Rind may be passed to `g`, `d`, `wx` and `r`, respectively.

see also: [`Grain`](@ref), [`Detrital`](@ref), [`WxAuth`](@ref), [`Rind`](@ref)

"""
function calcU(a::Vector,t::Vector;cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    cUout = zeros(Float64, length(a),length(t))
    cwouts = Vector{NamedTuple}(undef,length(a))
    for i in eachindex(a)
        isempty(cU) || (d.cU=cU[i])
        cwouts[i] = comminweath(a[i],g,d,wx,r)
        for j in eachindex(t)
            cUout[i,j] = drawdate(t[j],cwouts[i]).cU
        end
    end
    cUout
end