[![Build Status](https://travis-ci.org/flavioluiz/PFEM-article-supplementary-material.svg?branch=master)](https://travis-ci.org/flavioluiz/PFEM-article-supplementary-material)

# A Partitioned Finite Element Method for power-preserving discretization of open systems of conservation laws - Supplementary material
This archive contains supplementary material for the paper "A Partitioned Finite Element Method for power-preserving discretization of open systems of conservation laws", containing all the source codes for the numerial results presented in the paper.**  

The following codes are provided:
`codes/simulation1D_small.jl`: small amplitudes (linear) 1D simulation
`codes/simulation1D_large.jl`: large amplitudes (nonlinear) 1D simulation
`codes/simulation2D.jl`: large amplitudes (nonlinear) 2D simulation

## Installing
Before running the codes, you should install the [PortHamiltonian.jl](https://github.com/flavioluiz/PortHamiltonian.jl) package for Julia (it should work on v0.7):

```julia
julia> Pkg.clone("git://github.com/flavioluiz/PortHamiltonian.jl")
```

## 1D simulations

Two codes are provided for 1D simulations using PFEM, one using small amplitudes of boundary excitation, and the other with large amplitudes. It suffices to run the Julia file, by using:
```julia
julia> include("simulation1D_large.jl")
```
In the case of large amplitudes, the following results should be obtained (snapshots of the height at different times):

![Nonlinear simulations snapshots](./codes/simulations_snapshot_large.jpg)

![Nonlinear simulations time](./codes/fluid1DsimulationLarge.jpg)


## 2D simulations
One code is provided for 1D simulations using PFEM (`simulation2D.jl`). Changing the amplitude of the oscillation is straightforward (it suffices to change the variable "amp").

Example of results obtained:
