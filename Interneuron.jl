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

EL = -70*mV; t_i = 10*ms; C_i = 100*pF; v_thr = -50*mV

## Modelled as leaky-integrate-and-fire-neurons
dv_i_dt(v_i, I_ibg, W_EI, W_II, st_EI, st_II) = -(v_i .- EL) ./ t_i + (I_rec_i(W_EI, W_II, st_EI, st_II) .+ I_ibg) / C_i 

## External background current - uncorrelated activity 
## Constants 
mu = -100*pA; sigma = 400*pA; t_bg = 2*ms;

## Gaussian white noise with zero mean 
dI_ibg_dt(I_ibg) = -(I_ibg .- mu) ./ t_bg + sigma * rand(Normal(0.0, sqrt(dt)), size(I_ibg))

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons 
## Need update with individual types of interneurons 
## TODO - check that this is a correct calculation 
I_rec_i(W_EI, W_II, st_EI, st_II) = sum(abs.(W_EI) .* st_EI) .- sum(abs.(W_II) .* st_II)

## Implement spiking mechanism
## Update synaptic trace based on firing 
function simulateI(t, v_i, I_ibg, t_pyc, tspike_i, W_EI, W_II, st_EI, st_II)
    I_ibg += dI_ibg_dt(I_ibg) .* dt
    v_i += dv_i_dt(v_i, I_ibg, W_EI, W_II, st_EI, st_II) .* dt

    newt_ = map(x -> x >= v_thr ? 1 : 0, v_i)
    for i=1:length(tspike_i)
        if newt_[i] == 1 
            tspike_i[i] = t*dt
        end 
    end 
    map!(x -> x >= v_thr ? EL : x, v_i, v_i)

    return v_i, I_ibg, newt_, tspike_i
end

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end  
