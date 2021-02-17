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
dv_i_dt(v_i, I_ibg, u, R, W_EI, W_II, st_EI, st_II) = -(v_i .- EL) ./ t_i + (I_rec_i(u, R, W_EI, W_II, st_EI, st_II) + I_ibg) / C_i 

## External background current - uncorrelated activity 
## Constants 
u_i = -100*pA; theta_i = 400*pA; t_bg = 2*ms;

## Gaussian white noise with zero mean 
dI_ibg_dt(I_ibg) = add_gauss(-(I_ibg .- u_i) ./ t_bg, theta_i)

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons 
## Need update with individual types of interneurons 
I_rec_i(u, R, W_EI, W_II, st_EI, st_II) = sum.(abs.(W_EI) .* mu(u, R) .* st_EI) - sum.(abs.(W_II) .* st_II)

mu(u, R) = u .* R

## Constants
t_u = 100*ms; F = 0.1; t_R = 100*ms; 

## Define U 
du_dt(u, U) = -(u .- U) ./ t_u + (ones(size(u)) - u) .* F .* S
dR_dt(R, u) = -(R - ones(size(R))) ./ t_R - u .* R .* S

## Implement spiking mechanism
## Update synaptic trace based on firing 
function simulateI(t, v_i, I_ibg, t_, u, R, U, W_EI, W_II, st_EI, st_II, st_IE, st_II2)
    I_ibg += dI_ibg_dt(I_ibg) .* dt; #map!( x -> x < 0 ? 0 : x, I_ibg, I_ibg)
    v_i += dv_i_dt(v_i, I_ibg, u, R, W_EI, W_II, st_EI, st_II) .* dt

    u += du_dt(u, U) .* dt
    R += dR_dt(u, R) .* dt

    for i = 1:length(t_)
        if v_i[i] >= v_thr 
            # Update last spike timing 
            t_[i] = t*dt
            v_i[i] = EL
        end 

        # Update synaptic trace 
        st_IE[i] = update_st(t, t_[i], st_IE[i])
        st_II2[i] = update_st(t, t_[i], st_II2[i])
    end 

    return v_i, I_ibg, u, R, t_, st_IE, st_II2
end

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end  
