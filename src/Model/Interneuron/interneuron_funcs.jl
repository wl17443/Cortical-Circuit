#=============================================================================================================

Functions used for simulateI function

=============================================================================================================#

using Random, Distributions

## Background noise
dIbg_dt(Ibg, dt, noise_lvl) = - 0.5 .* (Ibg .+ 100) + noise_lvl["inh"] .* 400 .* rand(Normal(0.0, sqrt(dt)), size(Ibg))

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons
I_rec(st_EI, st_II, W_EI, W_II) = dropdims(sum(W_EI .* st_EI, dims=1), dims=1) - dropdims(sum(W_II .* st_II, dims=1), dims=1)
I_rec(st_EI, W_EI) = dropdims(sum(W_EI .* st_EI, dims=1), dims=1)

## Vip (IB) Interneurons
dv_vip_dt(v, u, I_bg, st_EI, st_II, W_EI, W_II) = (1.2 .* (v .+ 75) .* (v .+ 45) - u + I_rec(st_EI, st_II, W_EI, W_II) + I_bg) ./ 150
du_vip_dt(v, u) = 0.01 .* (5 .* (v .+ 75) - u)

## SST (LTS) Interneurons
dv_sst_dt(v, u, I_bg, st_EI, st_II, W_EI, W_II) = ((v .+ 56) .* (v .+ 42) - u + I_rec(st_EI, st_II, W_EI, W_II) + I_bg) ./ 100
du_sst_dt(v, u) = 0.03 .* (8 .* (v .+ 56) - u)

## PV (FS) Interneurons
dv_pv_dt(v, u, I_bg, st_EI, W_EI) = ((v .+ 55) .* (v .+ 40) - u + I_rec(st_EI, W_EI) + I_bg) ./ 20
du_pv_dt(v, u) = 0.2 .* (U.(v) - u)
U(v) = v < -55 ? 0 : 0.025 * (v + 55)^3
