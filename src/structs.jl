## struct support


"""
```julia
Detrital(cU, r48, 408)
```
`mutable struct` for initial detrital grain U-series parameters.
    
Constructor function `Detrital()` builds an instance with default values.

---
    
| field | units | description | 
| :----- | :-----: | :----------- |
| `cU`  | ng/g  | concentration of U in detrital phase | 
| `r48` |       | initial (²³⁴U/²³⁸U) of detrital grain |
| `r08` |       | initial (²³⁰Th/²³⁸U) of detrital grain | 

"""
mutable struct Detrital
    cU::Float64
    r48::Float64
    r08::Float64
end
Detrital() = Detrital(200,1,1)


"""

```julia
Grain(K, rfr, rho, L234, L230)
```
`mutable struct` for grain morphology and α-recoil parameters.
 
Constructor function `Grain()` builds an instance with default values.

---

| field | units | description | 
| :----- | :-----: | :----------- |
| `K`   |       | grain shape factorⁱ |
| `rfr` |       | surface roughness factorⁱⁱ |
| `rho` | kg/m³ | bulk density |
| `L234`| nm    | recoil length of ²³⁴U |
| `L230`| nm    | recoil length of ²³⁰Th |

---
ⁱ after Cartwright, 1962, doi:10.1093/annhyg/5.3.163, 
ⁱⁱ after White et al., 1996, doi: 10.1016/0016-7037(96)00106-8

"""
mutable struct Grain
    K::Float64
    rfr::Float64
    rho::Float64
    L234::Float64
    L230::Float64
end
Grain() = Grain(10,7,2650,34,37)

"""

```julia
Rind(evolve, p, z, rho, cU, r48, 408)
```
`mutable struct` for authigenic soluble rind parameters.
 
Constructor function `Rind()` builds an instance with default values.

---

| field | units | description | 
| :---- | :---: | :---------- |
| `evolve` | (boolean) | rind composition evolves (ON = `true`) |
| `p`   |       | proportion of grain covered by rind |
| `z`   |  μm   | thickness of rind |
| `rho` | kg/m³ | bulk density of rind phase |
| `cU`  | ng/g  | concentration of U in rind phase | 
| `r48` |       | initial (²³⁴U/²³⁸U) of rind phase|
| `r08` |       | initial (²³⁰Th/²³⁸U) of rind phase| 

"""
mutable struct Rind
    evolve::Bool
    p::Float64
    z::Float64
    rho::Float64
    cU::Float64
    r48::Float64
    r08::Float64
end
Rind() = Rind(true,0.4,0.1,2500,2000,3,4.6)

"""

```julia
WxAuth(sa_dependent, k, k_power, cU, r48, 408)
```
`mutable struct` for parameters of chemical weathering / authigenic replacement.
    
Constructor function `WxAuth()` builds an instance with default values.

---
    
| field | units | description | 
| :----- | :-----: | :----------- |
| `sa_dependent` | (boolean) | surface area dependence (ON = `true`)  |
| `k`   | g m⁻² a⁻¹ | rate of detrital phase alteration to authigenic phaseⁱ |
| `k_power` |   | exponent on `k` for (time-dependent) power-law weathering rate |
| `cU`  | ng/g  | concentration of U in authigenic phase | 
| `r48` |       | initial (²³⁴U/²³⁸U) of authigenic phase |
| `r08` |       | initial (²³⁰Th/²³⁸U) of authigenic phase | 
    
---
ⁱg m⁻² a⁻¹ if `sa_exponent`=`true`, (g g⁻¹) a⁻¹ if `sa_exponent`=`false`. This is an initial rate if `k_power` ≠ 1. 
ⁱⁱNot explored in this study, so it is set to 1 for a constant weathering rate, but it's there if you want it!

"""
mutable struct WxAuth
    sa_dependent::Bool
    k::Float64
    k_power::Float64
    cU::Float64
    r48::Float64
    r08::Float64
end
WxAuth() = WxAuth(true,1e-8,1,1000,3,4.6)