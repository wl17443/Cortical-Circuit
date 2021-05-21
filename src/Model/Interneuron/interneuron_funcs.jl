#=============================================================================================================

Functions used for simulateI function. 

=============================================================================================================#

using Random, Distributions

## Recurrent Inhibitory Inputs from Interneurons
## Interneurons get inhibitory input from other interneurons and excitatory input from connected pyramidal neurons
I_rec(st_E1I, st_E2I, W_EI) = dropdims(sum(W_EI .* st_E1I, dims=1), dims=1) + dropdims(sum(W_EI .* st_E2I, dims=1), dims=1)
I_rec(st_E1I, st_E2I, st_II, W_EI, W_II) = dropdims(sum(W_EI .* st_E1I, dims=1), dims=1) + dropdims(sum(W_EI .* st_E2I, dims=1), dims=1) - dropdims(sum(W_II .* st_II, dims=1), dims=1)

## Vip (IB) Interneurons
dv_vip_dt(v, u, st_E1I, st_E2I, st_II, W_EI, W_II, D, dt) = (1.2 .* (v .+ 75) .* (v .+ 45) - u + I_rec(st_E1I, st_E2I, st_II, W_EI, W_II) + D .* rand(Normal(0.0, sqrt(dt*1e-3)), length(v))) ./ 150
du_vip_dt(v, u) = 0.01 .* (5 .* (v .+ 75) - u)

## SST (LTS) Interneurons
dv_sst_dt(v, u, st_E1I, st_E2I, st_II, W_EI, W_II, D, dt) = ((v .+ 56) .* (v .+ 42) - u + I_rec(st_E1I, st_E2I, st_II, W_EI, W_II) + D .* rand(Normal(0.0, sqrt(dt*1e-3)), length(v))) ./ 100
du_sst_dt(v, u) = 0.03 .* (8 .* (v .+ 56) - u)

## PV (FS) Interneurons
dv_pv_dt(v, u, st_E1I, st_E2I, W_EI, D, dt) = ((v .+ 55) .* (v .+ 40) - u + I_rec(st_E1I, st_E2I, W_EI) + D .* rand(Normal(0.0, sqrt(dt*1e-3)), length(v))) ./ 20
du_pv_dt(v, u) = 0.2 .* (U.(v) - u)
U(v) = v < -55 ? 0 : 0.025 * (v + 55)^3
