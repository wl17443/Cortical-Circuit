## Pyramidal Neuron Model (3-Compartmental)

module PyramidalNeuron

include("Units.jl")
include("ModellingParameters.jl")
include("SomaticCompartment.jl")
include("DendriticCompartment.jl")
# include("AxonInitialSegment.jl")

using .Units
using .SomaticCompartment
using .DendriticCompartment
# using .AxonInitialSegment
using .ModellingParameters

v_thr = -50*mV

function simulatePyC(t, v_d, w_d, v_s, w_s, I_dbg, I_sbg, t_, st_SSTEd, st_PVEs)
    I_dbgPrime = I_dbg + dI_dbg_dt(I_dbg)*dt
    I_sbgPrime = I_sbg + dI_sbg_dt(I_sbg)*dt
    
    v_dPrime = v_d + dv_d_dt(v_d, I_dbg, w_d, last.(t_), st_SSTEd, t)*dt
    w_dPrime = w_d + dw_d_dt(w_d)*dt

    v_sPrime = v_s + dv_s_dt(v_s, I_sbg, w_s, st_PVEs)*dt
    w_sPrime = w_s + dw_s_dt(w_s, t)*dt

    ## TODO - If the new voltage surpasses the threshold, append to the spiking time and reset the voltage of the soma
    for i = 1:nr_pyc
        if v_sPrime[i] >= v_thr 
            append!(t_[i], t*dt)
            v_sPrime[i] = -70*mV
        end 
    end 

    return v_dPrime, w_dPrime, v_sPrime, w_sPrime, I_dbgPrime, I_sbgPrime, t_
end 

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 