## Pyramidal Neuron Model (3-Compartmental)

module PyramidalNeuron

include("Units.jl")
include("ModellingParameters.jl")
include("SomaticCompartment.jl")
include("DendriticCompartment.jl")
# include("AxonInitialSegment.jl")
include("UpdateSynapticTrace.jl")

using .Units
using .SomaticCompartment
using .DendriticCompartment
# using .AxonInitialSegment
using .ModellingParameters
using .UpdateSynapticTrace

v_thr = -50*mV

function simulatePyC(t, v_d, w_d, v_s, w_s, I_inj_d, I_inj_s, I_dbg, I_sbg, t_, W_SSTEd, W_PVEs, st_SSTEd, st_PVEs, st_EsSST, st_EsPV)
    I_dbg += dI_dbg_dt(I_dbg) .* dt; #map!( x -> x < 0 ? 0 : x, I_dbg, I_dbg)
    I_sbg += dI_sbg_dt(I_sbg) .* dt; #map!( x -> x < 0 ? 0 : x, I_sbg, I_sbg)
    
    v_d += dv_d_dt(v_d, I_inj_d, I_dbg, w_d, t_, W_SSTEd, st_SSTEd, t) .* dt
    w_d += dw_d_dt(w_d, v_d) .* dt

    v_s += dv_s_dt(v_s, v_d, I_inj_s, I_sbg, w_s, W_PVEs, st_PVEs) .* dt
    w_s += dw_s_dt(w_s, t) .* dt

    ## If the new voltage surpasses the threshold, append to the spiking time and reset the voltage of the soma
    for i = 1:nr_pyc
        if v_s[i] >= v_thr 
            t_[i] = t * dt
            v_s[i] = -70*mV
        end 

        st_EsSST[i] = update_st(t, t_[i], st_EsSST[i])
        st_EsPV[i] = update_st(t, t_[i], st_EsPV[i])
    end 

    return v_d, w_d, v_s, w_s, I_dbg, I_sbg, t_, st_EsSST, st_EsPV
end 

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 