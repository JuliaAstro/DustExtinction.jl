# DustExtinction.jl

```@contents
```

## Installation

From the REPL, press `]` to enter Pkg mode
```
(v 1.1) pkg> add DustExtinction
[...]

julia> using DustExtinction
```

## Usage

```jldoctest setup
julia> using DustExtinction
```

### Color laws

```jldoctest setup
julia> ccm89(4000., 3.1)
1.4645557029425842
```

These laws can be applied across higher dimension arrays using the `.` operator

```jldoctest setup
julia> ccm89.([4000., 5000.], 3.1)
2-element Array{Float64,1}:
 1.4645557029425842
 1.122246878899302

```

If you want to apply total extinction $A_V$ it's as simple as multiplcation
```jldoctest setup
julia> a_v=0.3
0.3

julia> a_v * ccm89(4000., 3.1)
0.43936671088277524
```


### Dust maps

```julia
julia> ENV["SFD98_DIR"] = "/home/user/data/dust"

# download maps (once)
julia> download_sfd98()

julia> dustmap = SFD98Map()
SFD98Map("/home/user/data/dust")

julia> ebv_galactic(dustmap, 0.1, 0.1)
0.793093095733043

julia> ebv_galactic(dustmap, [0.1, 0.2], [0.1, 0.2])
2-element Array{Float64,1}:
 0.793093
 0.539507
```


## Reference/API

### Extinction Laws
```@docs
ccm89
od94
cal00
```

### Dust Maps

```@docs
download_sfd98
SFD98Map
ebv_galactic
```
