module Connectivity

include("Units.jl")
include("ModellingParameters.jl")

using .Units
using .ModellingParameters
using Random, Distributions

## TODO - Weights are optimised with an optimisation method i.e. drift degree 

## Connectivity matrix that makes up the ring attractor network of PyC 
# W_EaEs = [ [0, 1, -1, -1, 1],
#            [1, 0, 1, -1, -1],
#            [-1, 1, 0, 1, -1],
#            [-1, -1, 1, 0, 1],
#            [1, -1, -1, 1, 0] ]
W_ESST = randn((nr_pyc, nr_sst))/(nr_pyc*nr_sst)
W_EPV = randn((nr_pyc, nr_pv))/(nr_pyc*nr_pv)

W_SSTEd = randn((nr_sst, nr_pyc))*0.2/(nr_sst*nr_pyc)
W_PVEs = randn((nr_pv, nr_pyc))*0.2/(nr_pv*nr_pyc)
# W_ChCEa = randn((nr_chc, nr_pyc))/(nr_chc*nr_pyc)

## Inter-Interneuron connectivity 
W_SSTPV = randn((nr_sst, nr_pv))/(nr_sst*nr_pv)
# W_SSTChC = randn((nr_sst, nr_chc))/(nr_sst*nr_chc)
# W_ChCPV = randn((nr_chc, nr_pv))/(nr_chc*nr_pv)
# W_ChCSST = randn((nr_chc, nr_sst))/(nr_chc*nr_sst)
W_PVSST = randn((nr_pv, nr_sst))/(nr_pv*nr_sst)
# W_PVChC = randn((nr_pv, nr_chc))/(nr_pv*nr_chc)

## U - Initial release probability 
U_ESST = randn(Uniform(0.1, 0.25), (nr_pyc, nr_sst))
U_EPV = randn(Uniform(0.1, 0.25), (nr_pyc, nr_pv))

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 