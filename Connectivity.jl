module Connectivity

include("Units.jl")
include("ModellingParameters.jl")

using .Units
using .ModellingParameters
using Random 

## Connectivity matrix that makes up the ring attractor network of PyC 
W_EaEs = [ [0, 1, -1, -1, 1],
           [1, 0, 1, -1, -1],
           [-1, 1, 0, 1, -1],
           [-1, -1, 1, 0, 1],
           [1, -1, -1, 1, 0] ]

W_SSTEd = randn((nr_sst, nr_pyc))/(nr_sst*nr_pyc)
W_PVEs = randn((nr_pv, nr_pyc))/(nr_pv*nr_pyc)
W_ChCEa = randn((nr_chc, nr_pyc))/(nr_chc*nr_pyc)

## Inter-Interneuron connectivity 
W_SSTPV = randn((nr_sst, nr_pv))/(nr_sst*nr_pv)
W_SSTChC = randn((nr_sst, nr_chc))/(nr_sst*nr_chc)
W_ChCPV = randn((nr_chc, nr_pv))/(nr_chc*nr_pv)
W_ChCSST = randn((nr_chc, nr_sst))/(nr_chc*nr_sst)
W_PVSST = randn((nr_pv, nr_sst))/(nr_pv*nr_sst)
W_PVChC = randn((nr_pv, nr_chc))/(nr_pv*nr_chc)

# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 