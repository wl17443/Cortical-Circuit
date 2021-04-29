using Random, Distributions

## Non-linear activation of the dendrite
E_d = -38e-3; D_d = 6e-3
f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants
EL = -70e-3; t_s = 16e-3; g_s = 1300e-9; C_s = 370e-9
b_s = -200e-9; t_s_w = 100e-3

## S(t) spike train from soma
dv_s_dt(v_s, v_d, I_inj, I_sbg, w_s, st_PVE, st_EE, W_PVE, W_EE) = -(v_s .- EL) ./ t_s + (g_s .* f(v_d) + w_s + I_inj + I_sbg + I_s_pv(st_PVE, W_PVE) + I_s_e(st_EE, W_EE)) ./ C_s
dw_s_dt(w_s, spike) = - w_s ./ t_s_w + b_s .* spike

## External background current - uncorrelated activity
## Constants
mu = 400e-9; t_bg = 2e-3; sigma = 450e-9;

## Gaussian white noise with zero mean and correlation
dI_sbg_dt(I_sbg, dt, noise_lvl) = -(I_sbg .- mu) ./ t_bg + noise_lvl["som"] .* sigma .* rand(Normal(0.0, sqrt(dt)), size(I_sbg))

## External input currents
I_s_pv(st_PVE, W_PVE) = -dropdims(sum(W_PVE .* st_PVE, dims=1), dims=1)
I_s_e(st_EE, W_EE) = dropdims(sum(W_EE .* st_EE, dims=1), dims=1)
