## Simple Model with PyC and INs
module RingAttractorModel

include("../model_components/PyramidalNeuron.jl")
include("../model_components/Interneuron.jl")

using ..Units
using .PyramidalNeuron: simulatePyc
using .Interneuron: simulateI
using Random, Distributions, Plots

export start_simulation

tau_syn = 5*ms

function start_simulation(t, dt, nr_pyc, nr_pv, nr_sst, bgnoise_lvl, alpha_exc, alpha_inh, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

steps = Int(t/dt)

## Variables
## Voltages
v_d = zeros(nr_pyc, steps); v_s = zeros(nr_pyc, steps)
w_d = zeros(nr_pyc, steps); w_s = zeros(nr_pyc, steps)
v_sst = zeros(nr_sst, steps); v_pv = zeros(nr_pv, steps)

## External background currents
I_sbg = zeros(nr_pyc, steps); I_dbg = zeros(nr_pyc, steps)
I_sstbg = zeros(nr_sst, steps); I_pvbg = zeros(nr_pv, steps)

## Spike train
t_pyc = zeros(nr_pyc, steps); t_sst = zeros(nr_sst, steps)
t_pv = zeros(nr_pv, steps)

## Last spike timing
t_soma = zeros(nr_pyc); tspike_sst = zeros(nr_sst)
tspike_pv = zeros(nr_sst)

## Synaptic trace
st_ESST = zeros(nr_pyc, nr_sst); st_EPV = zeros(nr_pyc, nr_pv)
st_SSTE = zeros(nr_sst, nr_pyc); st_PVE = zeros(nr_pv, nr_pyc)
st_SSTPV = zeros(nr_sst, nr_pv); st_PVSST = zeros(nr_pv, nr_sst)
st_EE = zeros(nr_pyc, nr_pyc)

## Initial values
## Voltages
v_d[1] = -70*mV; v_s[1] = -70*mV
v_sst[1] = -70*mV; v_pv[1] = -70*mV

## Simulation
for t = 2:steps
    v_d[:, t], w_d[:, t], v_s[:, t], w_s[:, t], I_dbg[:, t], I_sbg[:, t], t_pyc[:, t], t_soma[:] = simulatePyC(t, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], I_inj_d[:, t-1], I_inj_s[:, t-1], I_dbg[:, t-1], I_sbg[:, t-1], t_pyc[:, t-1], t_soma, st_SSTE, st_PVE, st_EE, W_SSTE, W_PVE, W_EE)

    ## Simulate for each type of interneuron
    v_sst[:, t], I_sstbg[:, t], t_sst[:, t], tspike_sst[:] = simulateI(t, v_sst[:, t-1], I_sstbg[:, t-1], tspike_sst, st_ESST, st_PVSST, W_ESST, W_PVSST)
    v_pv[:, t], I_pvbg[:, t], t_pv[:, t], tspike_pv[:] = simulateI(t, v_pv[:, t-1], I_pvbg[:, t-1], tspike_pv, st_EPV, st_SSTPV, W_EPV, W_SSTPV)

    ## Update synaptic trace
    for i=1:nr_pyc, j=1:nr_sst
        st_ESST[i, j] += (- st_ESST[i, j] ./ tau_syn) * dt + t_pyc[i, t]
    end
    for i=1:nr_pyc, j=1:nr_pv
        st_EPV[i, j] += (- st_EPV[i, j] ./ tau_syn) * dt + t_pyc[i, t]
    end
    for i=1:nr_sst, j=1:nr_pyc
        st_SSTE[i, j] += (- st_SSTE[i, j] ./ tau_syn) * dt + t_sst[i, t]
    end
    for i=1:nr_pv, j=1:nr_pyc
        st_PVE[i, j] += (- st_PVE[i, j] ./ tau_syn) * dt + t_pv[i, t]
    end
    for i=1:nr_sst, j=1:nr_pv
        st_SSTPV[i, j] += (- st_SSTPV[i, j] ./ tau_syn) * dt + t_sst[i, t]
    end
    for i=1:nr_pv, j=1:nr_sst
        st_PVSST[i, j] += (- st_PVSST[i, j] ./ tau_syn) * dt + t_pv[i, t]
    end
    for i=1:nr_pyc, j=1:nr_pyc
        st_EE[i, j] += (- st_EE[i, j] ./ tau_syn) * dt + t_pyc[i, t]
    end
end

## Plot results
step_list = [1:steps;]
display(Plots.heatmap(v_s, clim=(-70*mV, -50*mV), c=:balance))

end
end
