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
v_rest = -70*mV

function simulatePyC(t, v_d, w_d, v_s, w_s, I_inj_d, I_inj_s, I_dbg, I_sbg, t_, t_soma, W_SSTEd, W_PVEs, st_SSTEd, st_PVEs)
    I_dbg += dI_dbg_dt(I_dbg) .* dt
    I_sbg += dI_sbg_dt(I_sbg) .* dt
    
    v_d += dv_d_dt(v_d, I_inj_d, I_dbg, w_d, W_SSTEd, st_SSTEd, t_soma, t) .* dt
    w_d += dw_d_dt(w_d, v_d) .* dt

    v_s += dv_s_dt(v_s, v_d, I_inj_s, I_sbg, w_s, W_PVEs, st_PVEs) .* dt
    w_s += dw_s_dt(w_s, t_) .* dt

    newt_ = map(x -> x >= v_thr ? 1 : 0, v_s)
    for i=1:nr_pyc
        if newt_[i] == 1
            t_soma[i] = t*dt
        end 
    end 
    map!(x -> x >= v_thr ? v_rest : x, v_s, v_s)

    return v_d, w_d, v_s, w_s, I_dbg, I_sbg, newt_, t_soma 
end 

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 