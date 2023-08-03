## Measured Data

"""

```julia
ImpactChron.Tanaka2015()
```

Returns a NamedTuple of U-series data explicitly illustrating the effects of α-recoil implantation (Tanaka+ 2015, Chemical Geology, https://doi.org/10.1016/j.chemgeo.2014.12.025).

The data are from KOL48 (bulk) and magnetic susceptibility separates KOL48A–D (Table 3).

| fields | description |
| :----- | :---------- |
| r48 | (²³⁴U/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁴U/²³⁸U) | 
| r48 | (²³⁰Th/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁰Th/²³⁸U) | 

"""
Tanaka2015() = (
    r48 = [1.036, 1.049, 1.039, 1.032, 1.04],
    u48 = [0.0005, 0.001, 0.001, 0.0005, 0.001],
    r08 = [1.0425, 1.0586, 1.0392, 1.0350, 1.0642],
    u08 = [0.0027, 0.0054, 0.0039, 0.0021, 0.0022]
    )


"""

```julia
TaylorI()
```
Returns a NamedTuple of U-series data from detrital fractions of Taylor I drift sediments (measured in this study).

| fields | description |
| :----- | :---------- |
| r48 | (²³⁴U/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁴U/²³⁸U) | 
| r48 | (²³⁰Th/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁰Th/²³⁸U) | 

"""
TaylorI() = (
    r48 = [0.944502944, 0.955858312, 0.962708368, 0.951137072, 0.954462594, 0.964302593],
    u48 = [0.00116696, 0.002538465, 0.009153537, 0.006163414, 0.003866654, 0.001966957],
    r08= [0.943605822, 0.985127172, 1.068668456, 0.902913162, 1.026755662, 1.060088078],
    u08 = [0.020772331, 0.0223118, 0.025183063,	0.034591493, 0.041561336, 0.046794331] 
    )

"""

```julia
TaylorIII()
```
Returns a NamedTuple of U-series data from detrital fractions of Taylor III drift sediments (measured in this study).

| fields | description |
| :----- | :---------- |
| r48 | (²³⁴U/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁴U/²³⁸U) | 
| r48 | (²³⁰Th/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁰Th/²³⁸U) | 

"""
TaylorIII() = (
    r48 = [1.191904876, 1.200986111, 1.114035938, 1.184917544, 1.178136568, 1.091850317, 1.073923529, 1.105713217, 1.241336131],
    u48 = [0.010816397, 0.005698508, 0.034205831, 0.020006526, 0.01411507, 0.013214556, 0.029916289, 0.013160765, 0.016767997],
    r08 = [1.379008345, 1.365243125, 1.193035605, 1.244665406, 1.221846605, 1.160418328, 1.12382962, 1.125350103, 1.280323288],
    u08 = [0.026936737, 0.030721657, 0.025974771, 0.034491033, 0.034612668, 0.02803627, 0.030373848, 0.03289398, 0.03923534]
    )


"""

```julia
TaylorIV()
```
Returns a NamedTuple of U-series data from detrital fractions of Taylor IV drift sediments (measured in this study).

| fields | description |
| :----- | :---------- |
| r48 | (²³⁴U/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁴U/²³⁸U) | 
| r48 | (²³⁰Th/²³⁸U) |
| u48 | 2σ uncertainty of (²³⁰Th/²³⁸U) | 

"""
TaylorIV() = (
    r48 = [0.989544381, 0.988045986, 1.034253721, 1.043602168, 1.015417282, 1.112719734, 1.087193038, 1.106630613], 
    u48 = [0.001069865, 0.004102872, 0.003387656, 0.003343425, 0.002353557, 0.002495865, 0.007035458, 0.027496583], 
    r08 = [1.04123385, 1.005395862, 1.121889924, 1.125412136, 1.077348361, 1.256509829, 1.194628943, 1.228268553], 
    u08 = [0.02257361, 0.023122304, 0.037680406, 0.025937369, 0.026106258, 0.0520976, 0.047257529, 0.031318882] 
    )