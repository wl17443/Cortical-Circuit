#=============================================================================================================

Function to simulate PyC dendritic compartment

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
dv_d_dt(v, Iinj, Ibg, w_d, st_SSTE, W_SSTE, t_soma, t) = -(v .- EL) ./ t_d + (g_d .* f(v) + c_d .* K(t .- t_soma) + w_d + Iinj + Ibg + I_d_sst(st_SSTE, W_SSTE)) ./ C_d
dw_d_dt(w_d, v) = - w_d ./ t_d_w + a_d .* (v .- EL) ./ t_d_w

## External background current - uncorrelated activity
## Constants
mu = -300e-9; tau = 2e-3
sigma = 450e-9;

## Gaussian white noise with zero mean with variance dt
dI_dbg_dt(Ibg, dt, noise_lvl) = -(Ibg .- mu) ./ tau + noise_lvl["den"] .* sigma .* rand(Normal(0.0, sqrt(dt)), size(Ibg))

## External input current
I_d_sst(st_SSTE, W_SSTE) = -dropdims(sum(W_SSTE .* st_SSTE, dims=1), dims=1)
