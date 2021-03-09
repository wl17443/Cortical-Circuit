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

# W_EE = rand(Normal(0.37, 0.1*sqrt(16)), (nr_pyc, nr_pyc)) .* nA
W_EE = zeros(nr_pyc, nr_pyc)
for i=1:nr_pyc
    W_EE[i, i+1 == nr_pyc+1 ? 1 : i+1] = abs(rand(Normal(0.37, 0.1*sqrt(16)))) * nA
    W_EE[i, i-1 == 0 ? nr_pyc : i-1] = abs(rand(Normal(0.37, 0.1*sqrt(16)))) * nA
end

## Synaptic Weights (Fixed)
## E->I{SST,PV}
W_ESST = rand(Normal(0.37, 0.1*sqrt(11)), (nr_pyc, nr_sst)) .* nA
W_EPV = rand(Normal(0.82, 0.1*sqrt(23)), (nr_pyc, nr_pv)) .* nA

## I{SST,PV}->E{s,d}
W_SSTEd = rand(Normal(0.49, 0.11*sqrt(20)), (nr_sst, nr_pyc)) .* nA
W_PVEs = rand(Normal(0.52, 0.11*sqrt(21)), (nr_pv, nr_pyc)) .* nA

## I{SST,PV}->I{SST,PV}
W_SSTPV = rand(Normal(0.37, 0.1*sqrt(11)), (nr_sst, nr_pv)) .* nA
W_PVSST = rand(Normal(0.83, 0.25*sqrt(7)), (nr_pv, nr_sst)) .* nA

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end
