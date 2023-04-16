## Decay Constants, in a⁻¹

# ²³⁸U --- Jaffey et al. 1971, Phys. Rev. C, doi: 10.1103/PhysRevC.4.1889
const l238 = log(2)/ 4.4683e9

# ²³⁴U, ²³⁰Th --- Cheng et al. 2013, EPSL, doi: 10.1016/j.epsl.2013.04.006
const l234 = log(2)/ 245620 
const l230 = log(2)/ 75584

# ²³²Th --- Audi et al. 1997, Nuc. Phys. A, doi: 10.1016/S0375-9474(97)00482-X
const l232 = log(2)/ 14.05e9

# ²²⁶Ra --- Holden, 1990, Pure and Applied Chem., doi: 10.1351/pac199062050941
const l226 = log(2)/ 1599 

## Secular equilibrium values based on decay constants:
const seU = l238/l234
const seThU = l238/l230