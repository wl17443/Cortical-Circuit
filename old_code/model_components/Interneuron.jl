## Interneuron model (Leaky-Integrate-And-Fire)
module Interneuron

using ..Units
using Random, Distributions
using Noise

export simulateI

## Modelled as leaky-integrate-and-fire-neurons
EL = -70*mV; tau = 10*ms; C = 100*pF; v_thr = -65*mV
dv_dt(v, I_bg, st_EI, st_II, W_EI, W_II) = -(v .- EL) ./ tau + (I_rec(st_EI, st_II, W_EI, W_II) + I_bg) ./ C

## External background current - uncorrelated activity
## Constants
mu = -100*pA; t_bg = 2*ms;
sigma = 400*pA;

## Gaussian white noise with zero mean
dIbg_dt(Ibg) = -(Ibg .- mu) ./ t_bg + bgnoise_lvl .* sigma .* rand(Normal(0.0, sqrt(dt)), size(Ibg))

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons
## Need update with individual types of interneurons
I_rec(st_EI, st_II, W_EI, W_II) = sum(W_EI * st_EI, dims=2) - sum(W_II * st_II, dims=2)

## Implement spiking mechanism
## Update synaptic trace based on firing
function simulateI(t, v, Ibg, tspike, st_EI, st_II, W_EI, W_II)
    v += dv_dt(v, Ibg, st_EI, st_II, W_EI, W_II) .* dt
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

end
