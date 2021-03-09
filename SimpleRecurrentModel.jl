## Simple Model with PyC and INs
module SimpleRecurrentNetwork

include("PyramidalNeuron.jl")
include("Interneuron.jl")
include("Units.jl")
include("ModellingParameters.jl")
# include("Connectivity.jl")

using .PyramidalNeuron
using .Interneuron
using .Units
using .ModellingParameters
using .Connectivity
using Random, Distributions
using Plots
using Noise
using Dates

## Variables
## Voltages
v_d = zeros(nr_pyc, steps)
w_d = zeros(nr_pyc, steps)
v_s = zeros(nr_pyc, steps)
w_s = zeros(nr_pyc, steps)
v_sst = zeros(nr_sst, steps)
v_pv = zeros(nr_pv, steps)

## External background currents
I_sbg = zeros(nr_pyc, steps)
I_dbg = zeros(nr_pyc, steps)
I_sstbg = zeros(nr_sst, steps)
I_pvbg = zeros(nr_pv, steps)

## Spike train
t_pyc = zeros(nr_pyc, steps)
t_sst = zeros(nr_sst, steps)
t_pv = zeros(nr_pv, steps)

## Last spike timing
t_soma = zeros(nr_pyc)
tspike_sst = zeros(nr_sst)
tspike_pv = zeros(nr_sst)

## Synaptic trace
st_ESST = zeros(nr_pyc, nr_sst)
st_EPV = zeros(nr_pyc, nr_pv)
st_SSTE = zeros(nr_sst, nr_pyc)
st_PVE = zeros(nr_pv, nr_pyc)
st_SSTPV = zeros(nr_sst, nr_pv)
st_PVSST = zeros(nr_pv, nr_sst)
st_EE = zeros(nr_pyc, nr_pyc)
## Initial values
v_d[1] = -70*mV
v_s[1] = -70*mV
v_sst[1] = -70*mV
v_pv[1] = -70*mV

## Injected Current
I_inj_d = zeros(nr_pyc, steps)
I_inj_s = zeros(nr_pyc, steps)

## Setup input current
I_inj_s[3, 1:100] .= 10*nA
I_inj_s[4, 5000:5100] .= 10*nA

## Simulation
for t = 2:steps
    v_d[:, t], w_d[:, t], v_s[:, t], w_s[:, t], I_dbg[:, t], I_sbg[:, t], t_pyc[:, t], t_soma[:] = simulatePyC(t, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], I_inj_d[:, t-1], I_inj_s[:, t-1], I_dbg[:, t-1], I_sbg[:, t-1], t_pyc[:, t-1], t_soma, st_SSTE, st_PVE, st_EE)

    ## Simulate for each type of interneuron
    v_sst[:, t], I_sstbg[:, t], t_sst[:, t], tspike_sst[:] = simulateI("SST", t, v_sst[:, t-1], I_sstbg[:, t-1], tspike_sst, st_ESST, st_PVSST)
    v_pv[:, t], I_pvbg[:, t], t_pv[:, t], tspike_pv[:] = simulateI("PV", t, v_pv[:, t-1], I_pvbg[:, t-1], tspike_pv, st_EPV, st_SSTPV)

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

    # st_EsSST[:] += (- st_EsSST ./ tau_syn) .* dt + t_pyc[:, t]
    # map!(((x, i), ) -> map!(y -> y += (- y / tau_syn) * dt + t_pyc[i, t], x, x), st_EsSST, st_EsSST)
    # st_EsPV[:] += (- st_EsPV ./ tau_syn) .* dt + t_pyc[:, t]
    # st_SSTEd[:] += (- st_SSTEd ./ tau_syn) .* dt + t_sst[:, t]
    # st_PVEs[:] += (- st_PVEs ./ tau_syn) .* dt + t_pv[:, t]
    # st_SSTPV[:] += (- st_SSTPV ./ tau_syn) .* dt + t_sst[:, t]
    # st_PVSST[:] += (- st_PVSST ./ tau_syn) .* dt + t_pv[:, t]
end

step_list = [1:steps;]
# p1 = plot(step_list, v_s[:], label="Somatic Voltage")
# p2 = plot(step_list, v_d[:], label="Dendritic Voltage")
# p3 = plot(step_list, v_sst[:], label="SST Voltage")
# p4 = plot(step_list, v_pv[:], label="PV Voltage")

# p5 = plot(step_list, I_sbg[:], label="Somatic Background Current")
# p6 = plot(step_list, I_dbg[:], label="Dendritic Background Current")
# p7 = plot(step_list, I_sstbg[:], label="SST Background Current")
# p8 = plot(step_list, I_pvbg[:], label="PV Background Current")

# p9 = plot(step_list, t_pyc[:], label="PyC Spike Train")
# p10 = plot(step_list, t_sst[:], label="SST Spike Train")
# p11 = plot(step_list, t_pv[:], label="PV Spike Train")

p12 = plot(step_list, t_pyc[1, :], label="PyC1")
p13 = plot(step_list, t_pyc[2, :], label="PyC2")
p14 = plot(step_list, t_pyc[3, :], label="PyC3")
p15 = plot(step_list, t_pyc[4, :], label="PyC4")
p16 = plot(step_list, t_pyc[5, :], label="PyC5")

# Plot voltage trace
# plt = plot(p1, p3, p4, layout=(3,1))
# plt2 = plot(p9, p10, p11, layout=(3,1))
# datetime = now()
# fdatetime = Dates.format(datetime, "yyyy-mm-dd-HH-MM-SS")
# savefig(plt, "./fig/voltage-trace-$fdatetime.png")
# savefig(plt2, "./fig/spike-train-$fdatetime.png")

# Plot voltage trace and background current
# display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(4,2), size=(1000,500)))

# Plot spike train
# display(plot(p1, p2, p3, p4, p9, p10, p11))

# Plot voltage train of ring attractor
display(plot(p12, p13, p14, p15, p16, layout=(5,1)))
end
