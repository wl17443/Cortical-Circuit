#=============================================================================================================

Function to simulate PyC dendritic compartment. 

=============================================================================================================#

using Random, Distributions

## Non-linear activation of the dendrite
E_d = -38e-3; D_d = 6e-3
f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants
EL = -70e-3; t_d = 7e-3; g_d = 1200e-9; c_d = 2600e-9; t_d_w = 30e-3; a_d = -13e-6; C_d = 170e-9

## Boxcar kernel K
function K(t)
    for i = 1:length(t)
        if t[i] < 0.5e-3
            t[i] = 0
        elseif t[i] >= 0.5e-3 && t[i] < 2.5e-3
            t[i] = 1
        else t[i] = 0 end
    end
    return t
end

##  where t_ is the last spike time of soma - updates with every spike (global variable)
dv_d_dt(v, Iinj, w_d, st_IE, W_IE, t_soma, t, D, dt) = -(v .- EL) ./ t_d + (g_d .* f(v) + c_d .* K(t .- t_soma) + w_d + Iinj + Irec(st_IE, W_IE) + D .* rand(Normal(0.0, sqrt(dt)), length(v))) ./ C_d
dv_d_dt(v, w_d, st_EE, st_IE, W_EE, W_IE, t_soma, t, D, dt) = -(v .- EL) ./ t_d + (g_d .* f(v) + c_d .* K(t .- t_soma) + w_d + Irec(st_EE, st_IE, W_EE, W_IE) + D .* rand(Normal(0.0, sqrt(dt)), length(v))) ./ C_d
dw_d_dt(w_d, v) = - w_d ./ t_d_w + a_d .* (v .- EL) ./ t_d_w

## External input current
Irec(st_IE, W_IE) = - dropdims(sum(W_IE .* st_IE, dims=1), dims=1)
Irec(st_EE, st_IE, W_EE, W_IE) = dropdims(sum(W_EE .* st_EE, dims=1), dims=1) - dropdims(sum(W_IE .* st_IE, dims=1), dims=1)
