## Dendritic Compartment
module DendriticCompartment

include("Units.jl")
include("ModellingParameters.jl")
include("Connectivity.jl")

using .Units
using .ModellingParameters
using .Connectivity
using Random, Distributions 

## Non-linear activation of the dendrite
E_d = -38*mV; D_d = 6*mV 

f(v) = 1 ./ (1 .+ exp.(-(v .- E_d) ./ D_d))

## Constants 
EL = -70*mV; t_d = 7*ms; g_d = 1200*pA; c_d = 2600*pA; t_d_w = 30*ms; a_d = -13*nS; C_d = 170*pF

## Boxcar kernel K 
function K(t)
    for i = 1:size(t)[1]
        if t[i] < 1*ms
            t[i] = 0  
        elseif t[i] > 1*ms && t[i] < 3*ms
            t[i] = -50*mV ## TODO - Change to actual plateau value
        else t[i] = 0 end 
    end 
    return t
end 

##  where t_ is the last spike time of soma - updates with every spike (global variable)
## TODO - Check that f(v_d) is not f(v_s)
## TODO - Define I_d 
dv_d_dt(v_d, I_dbg, w_d, t_, st_SSTEd, t) = -(v_d .- EL) ./ t_d + (g_d .* f(v_d) + c_d .* K(t*dt .- t_) + w_d + I_dbg + I_d_sst(st_SSTEd)) ./ C_d
dw_d_dt(w_d, v_d) = - w_d ./ t_d_w + a_d .* (v_d .- EL) ./ t_d_w

## External background current - uncorrelated activity 
## Constants 
u_d = 400*pA; theta_d = 450*pA; t_bg = 2*ms

## Gaussian white noise with zero mean 
dI_dbg_dt(I_dbg) = -(I_dbg .- u_d) ./ t_bg + theta_d .* randn(size(I_dbg))

I_d_sst(st_SSTEd) = -sum.(abs.(W_SSTEd)*st_SSTEd)

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 