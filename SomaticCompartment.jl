## Somatic compartment
module SomaticCompartment

include("Units.jl")

using .Units 
using Random, Distributions 

## Non-linear activation of the dendrite
E_d = -38*mV; D_d = 6*mV 

f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants 
EL = -70*mV; t_s = 16*ms; g_s = 1300*pA; C_s = 370*pF
b_s = -300*pA; t_s_w = 100*ms 

## S(t) spike train from soma 
dv_s_dt(v_s, I_sbg, t) = -(v_s .- EL) ./ t_s + (g_s .* f(v_d) .+ w_s + I_sbg + I_s)/C_s
dw_s_dt(w_s, spike) = - w_s ./ t_s_w .+ b_s * spike

## External background current - uncorrelated activity 
## Constants 
u_s = 400*pA; theta_s = 450*pA; t_bg = 2*ms

## Gaussian white noise with zero mean and correlation
dI_sbg_dt(I_sbg, t) = -(I_sbg .- u_s) ./ t_bg + theta_s .* randn(size(I_sbg))

# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 