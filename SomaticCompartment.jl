## Somatic compartment
module SomaticCompartment

include("Units.jl")
include("Connectivity.jl")
include("ModellingParameters.jl")

using .Units 
using .Connectivity
using .ModellingParameters
using Random, Distributions 
using Noise 

## Non-linear activation of the dendrite
E_d = -38*mV; D_d = 6*mV 

f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants 
EL = -70*mV; t_s = 16*ms; g_s = 1300*pA; C_s = 370*pF
b_s = -200*pA; t_s_w = 100*ms 

## S(t) spike train from soma 
dv_s_dt(v_s, v_d, I_inj, I_sbg, w_s, st_PVEs) = -(v_s .- EL) ./ t_s + (g_s .* f(v_d) .+ w_s .+ I_inj .+ I_sbg .+ I_s_pv(st_PVEs)) ./ C_s
dw_s_dt(w_s, spike) = - w_s ./ t_s_w .+ b_s .* spike

## External background current - uncorrelated activity 
## Constants 
mu = 400*pA; sigma = 450*pA; t_bg = 2*ms

## Gaussian white noise with zero mean and correlation
dI_sbg_dt(I_sbg) = -(I_sbg .- mu) ./ t_bg + sigma * rand(Normal(0.0, sqrt(dt)), size(I_sbg))

I_s_pv(st_PVEs) = -sum(abs.(W_PVEs) * st_PVEs)

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 