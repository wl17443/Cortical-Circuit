## Interneuron model (Leaky-Integrate-And-Fire)

module Interneuron

include("Units.jl")
# include("Connectivity.jl")
include("ModellingParameters.jl")
include("UpdateSynapticTrace.jl")

using .Units
# using .Connectivity
using .ModellingParameters
using .UpdateSynapticTrace
using Random, Distributions 
using Noise 

EL = -70*mV; tau = 10*ms; C = 100*pF; v_thr = -50*mV

## Modelled as leaky-integrate-and-fire-neurons
dv_dt(v, I_bg, W_EI, W_II, st_EI, st_II) = -(v .- EL) ./ tau + (I_rec(W_EI, W_II, st_EI, st_II) .+ I_bg) / C

## External background current - uncorrelated activity 
## Constants 
mu = -100*pA; sigma = 400*pA; t_bg = 2*ms;

## Gaussian white noise with zero mean 
dIbg_dt(Ibg) = -(Ibg .- mu) ./ t_bg + sigma * rand(Normal(0.0, sqrt(dt)), size(Ibg))

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons 
## Need update with individual types of interneurons 
## TODO - check that this is a correct calculation 
I_rec(W_EI, W_II, st_EI, st_II) = sum(abs.(W_EI) .* st_EI) - sum(abs.(W_II) .* st_II)

## Implement spiking mechanism
## Update synaptic trace based on firing 
function simulateI(t, v, Ibg, t_pyc, tspike, W_EI, W_II, st_EI, st_II)
    v += dv_dt(v, Ibg, W_EI, W_II, st_EI, st_II) .* dt
    Ibg += dIbg_dt(Ibg) .* dt

    newt_ = map(x -> x >= v_thr ? 1 : 0, v)
    for i=1:length(tspike)
        if newt_[i] == 1 
            tspike[i] = t*dt
        end 
    end 
    map!(x -> x >= v_thr ? EL : x, v, v)

    return v, Ibg, newt_, tspike
end

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end  
