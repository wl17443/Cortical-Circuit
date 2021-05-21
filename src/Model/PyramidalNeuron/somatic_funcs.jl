#=============================================================================================================

Functions to simulate PyC stomatic compartment. 

=============================================================================================================#

using Random, Distributions

## Non-linear activation of the dendrite
E_d = -38e-3; D_d = 6e-3
f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants
EL = -70e-3; t_s = 16e-3; g_s = 1300e-9; C_s = 370e-9
b_s = -200e-9; t_s_w = 100e-3

## S(t) spike train from soma
# Layer 1
dv1_s_dt(v, v_d, w_s, st_IE, st_EE, st_E2E1, W_IE, W_EE, D, dt) = -(v .- EL) ./ t_s + (g_s .* f(v_d) + w_s + Irec(st_EE, st_E2E1, st_IE, W_EE, W_IE) + D .* rand(Normal(0.0, sqrt(dt)), length(v))) ./ C_s
# Layer 2
dv2_s_dt(v, v_d, I_inj, w_s, st_IE, st_EE, W_IE, W_EE, D, dt) = -(v .- EL) ./ t_s + (g_s .* f(v_d) + w_s + I_inj + Irec(st_EE, st_IE, W_EE, W_IE) + D .* rand(Normal(0.0, sqrt(dt)), length(v))) ./ C_s
dw_s_dt(w_s, spike) = - w_s ./ t_s_w + b_s .* spike

## External input currents
Irec(st_EE, st_IE, W_EE, W_IE) = dropdims(sum(W_EE .* st_EE, dims=1), dims=1) - dropdims(sum(W_IE .* st_IE, dims=1), dims=1)
Irec(st_EE, st_E2E1, st_IE, W_EE, W_IE) = dropdims(sum(W_EE .* st_EE, dims=1), dims=1) + dropdims(sum(W_EE .* st_E2E1, dims=1), dims=1) - dropdims(sum(W_IE .* st_IE, dims=1), dims=1)
