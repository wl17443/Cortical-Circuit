## Pyramidal Neuron Model (3-Compartmental)

module PyramidalNeuron

include("Units.jl")
include("SomaticCompartment.jl")
include("DendriticCompartment.jl")

using .Units
using .SomaticCompartment
using .DendriticCompartment

function simulatePyC(t_step, dt, v_d, w_d, v_s, w_s, I_dbg, I_sbg)
    I_dbgPrime = I_dbg + dI_dbg_dt(I_dbg, t_step)*dt
    I_sbgPrime = I_sbg + dI_sbg_dt(I_sbg, t_step)*dt
    
    v_dPrime = v_d + dv_d_dt(v_d, I_dbgPrime, t_step)*dt
    w_dPrime = w_d + dw_d_dt(w_d, t_step)*dt

    v_sPrime = v_s + dv_s_dt(v_s, I_sbgPrime, t_step)*dt
    w_sPrime = w_s + dw_s_dt(w_s, t_step)*dt

    return v_dPrime, w_dPrime, v_sPrime, w_dPrime 
end 

# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 