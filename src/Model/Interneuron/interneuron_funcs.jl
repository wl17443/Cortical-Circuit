#=============================================================================================================

Functions used for simulateI function

=============================================================================================================#

using Random, Distributions

## Modelled as leaky-integrate-and-fire-neurons
EL = -70e-3; tau = 10e-3; C = 100e-9; v_thr = -65e-3
dv_dt(v, I_bg, st_EI, st_II, W_EI, W_II) = -(v .- EL) ./ tau + (I_rec(st_EI, st_II, W_EI, W_II) + I_bg) ./ C

## External background current - uncorrelated activity
mu = -100e-9; t_bg = 2e-3; sigma = 400e-9;
dIbg_dt(Ibg, dt, noise_lvl) = -(Ibg .- mu) ./ t_bg + noise_lvl .* sigma .* rand(Normal(0.0, sqrt(dt)), size(Ibg))

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons
I_rec(st_EI, st_II, W_EI, W_II) = dropdims(sum(W_EI .* st_EI, dims=1), dims=1) - dropdims(sum(W_II .* st_II, dims=1), dims=1)
