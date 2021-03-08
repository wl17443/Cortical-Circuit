## Dendritic Compartment
module DendriticCompartment

include("Units.jl")
include("ModellingParameters.jl")
include("Connectivity.jl")

using .Units
using .ModellingParameters
using .Connectivity
using Random, Distributions 
using Noise 

## Non-linear activation of the dendrite
E_d = -38*mV; D_d = 6*mV 

f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants 
EL = -70*mV; t_d = 7*ms; g_d = 1200*pA; c_d = 2600*pA; t_d_w = 30*ms; a_d = -13*nS; C_d = 170*pF

## Boxcar kernel K 
function K(t)
    for i = 1:length(t)
        if t[i] < 0.5*ms
            t[i] = 0
        elseif t[i] >= 0.5*ms && t[i] < 2.5*ms
            t[i] = 1
        else t[i] = 0 end 
    end 
    return t
end 

##  where t_ is the last spike time of soma - updates with every spike (global variable)
dv_d_dt(v_d, I_inj, I_dbg, w_d, st_SSTEd, t_soma, t) = -(v_d .- EL) ./ t_d + (g_d .* f(v_d) .+ c_d .* K(t .- t_soma) .+ w_d .+ I_inj .+ I_dbg .+ I_d_sst(st_SSTEd)) ./ C_d
dw_d_dt(w_d, v_d) = - w_d ./ t_d_w + a_d .* (v_d .- EL) ./ t_d_w

## External background current - uncorrelated activity 
## Constants 
u_d = -300*pA; theta_d = 450*pA; t_bg = 2*ms

## Gaussian white noise with zero mean with variance dt 
dI_dbg_dt(I_dbg) = -(I_dbg .- u_d) ./ t_bg + rand(Normal(0.0, sqrt(dt)), size(I_dbg))

I_d_sst(st_SSTEd) = -sum(abs.(W_SSTEd) * st_SSTEd)

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 