## Pyramidal Neuron Model (3-Compartmental)

module PyramidalNeuron

include("Units.jl")
include("SomaticCompartment.jl")
include("DendriticCompartment.jl")
include("AxonInitialSegment.jl")

using .Units
using .SomaticCompartment
using .DendriticCompartment
using .AxonInitialSegment

function simulatePyC(t_step, dt, v_d, w_d, v_s, w_s, I_dbg, I_sbg, syntrace_sst, syntrace_pv)
    I_dbgPrime = I_dbg + dI_dbg_dt(I_dbg)*dt
    I_sbgPrime = I_sbg + dI_sbg_dt(I_sbg)*dt
    
    v_dPrime = v_d + dv_d_dt(v_d, I_dbgPrime, t_step)*dt
    w_dPrime = w_d + dw_d_dt(w_d)*dt

    v_sPrime = v_s + dv_s_dt(v_s, I_sbgPrime, syntrace_pv)*dt
    w_sPrime = w_s + dw_s_dt(w_s, t_step)*dt

    return v_dPrime, w_dPrime, v_sPrime, w_dPrime 
end 

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 