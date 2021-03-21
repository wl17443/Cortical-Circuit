## Somatic compartment
module SomaticCompartment

using ..Units
using Random, Distributions
using Noise

export dv_s_dt, dw_s_dt, dI_sbg_dt

## Non-linear activation of the dendrite
E_d = -38*mV; D_d = 6*mV
f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants
EL = -70*mV; t_s = 16*ms; g_s = 1300*pA; C_s = 370*pF
b_s = -200*pA; t_s_w = 100*ms

## S(t) spike train from soma
dv_s_dt(v_s, v_d, I_inj, I_sbg, w_s, st_PVE, st_EE, W_PVE, W_EE) = -(v_s .- EL) ./ t_s + (g_s .* f(v_d) + w_s + I_inj + I_sbg + I_s_pv(st_PVE, W_PVE) + I_s_e(st_EE, W_EE)) ./ C_s
dw_s_dt(w_s, spike) = - w_s ./ t_s_w + b_s .* spike

## External background current - uncorrelated activity
## Constants
mu = 400*pA; t_bg = 2*ms
sigma = 450*pA;

## Gaussian white noise with zero mean and correlation
dI_sbg_dt(I_sbg) = -(I_sbg .- mu) ./ t_bg + bgnoise_lvl .* sigma .* rand(Normal(0.0, sqrt(dt)), size(I_sbg))

## External input currents
I_s_pv(st_PVE, W_PVE) = -sum(W_PVE * st_PVE, dims=2)
I_s_e(st_EE, W_EE) = sum(W_EE * st_EE, dims=2)

end
