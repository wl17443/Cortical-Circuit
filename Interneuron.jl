## Interneuron model (Leaky-Integrate-And-Fire)

module Interneuron

include("Units.jl")
include("Connectivity.jl")

using .Units
using .Connectivity
using Random, Distributions 

EL = -70*mV; t_i = 10*ms; C_i = 100*pF;

## Modelled as leaky-integrate-and-fire-neurons
dv_i_dt(v_i, I_ibg, u, R, t) = -(v_i .- EL) ./ t_i + (I_rec_i(u, R, t) + I_ibg) / C_i 

## External background current - uncorrelated activity 
## Constants 
u_i = -100*pA; theta_i = 400*pA; t_bg = 2*ms;

## Gaussian white noise with zero mean 
dI_ibg_dt(I_ibg) = -(I_ibg .- u_i) ./ t_bg + theta_i .* randn(size(I_ibg))

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons 
## TODO - Need update with individual types of interneurons 
I_rec_i(u, R, t_step) = W_EI * u(u, R) * s(t_step) - W_II * s(t_step)

u(u, R) = u * R

## Constants
t_u = 100*ms; F = 0.1; S = 1; t_R = 100*ms; 

du_dt(u) = -(u .- U) ./ t_u + (ones(size(u)) - u) .* F .* S
dR_dt(R, u) = -(R - ones(size(R))) ./ t_R - u .* R .* S

function simulateI(t_step, dt, v_i, I_ibg, u, R)
    I_ibgPrime = I_ibg + dI_ibg_dt(I_ibg) * dt
    uPrime = u + du_dt(u) .* dt
    RPrime = R + dR_dt(u, R) .* dt
    v_iPrime = dv_i_dt(v_i, I_ibg, u, R, t_step)

    return V_iPrime, I_ibgPrime
end  

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end  
