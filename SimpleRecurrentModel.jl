## Simple Model with PyC and INs 
module SimpleRecurrentNetwork 

include("PyramidalNeuron.jl")
include("Interneuron.jl")
# include("Connectivity.jl")
include("Units.jl")
include("ModellingParameters.jl")

using .PyramidalNeuron
using .Interneuron
# using .Connectivity
using .Units    
using .ModellingParameters
using Random, Distributions
using Plots 
using Noise 

initial_synaptic_weight = 4*nS

## Weights being optimised 
## E->I{SST,PV}
W_ESST = zeros((nr_pyc, nr_sst)) 
W_EPV = zeros((nr_pyc, nr_pv)) 

## I{SST,PV}->E{s,d}
W_SSTEd = zeros((nr_sst, nr_pyc)) 
W_PVEs = zeros((nr_pv, nr_pyc)) 

## I{SST,PV}->I{SST,PV}
W_SSTPV = zeros((nr_sst, nr_pv))
W_PVSST = zeros((nr_pv, nr_sst))

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
st_EsSST = zeros(nr_pyc, nr_sst)
st_EsPV = zeros(nr_pyc, nr_pv)
st_SSTEd = zeros(nr_sst, nr_pyc)
st_PVEs = zeros(nr_pv, nr_pyc)
st_SSTPV = zeros(nr_sst, nr_pv)
st_PVSST = zeros(nr_pv, nr_sst)

## Initial values 
v_d[1] = -70*mV
v_s[1] = -70*mV
v_sst[1] = -70*mV
v_pv[1] = -70*mV

## Injected Current 
I_inj_d = zeros(nr_pyc, steps)
I_inj_s = zeros(nr_pyc, steps)

u_bar_minus = zeros(nr_pyc)
u_bar_plus = zeros(nr_pyc)

tau_u_bar_minus = 7*ms; tau_u_bar_plus = 10*ms 
## Where u is the postsynaptic membrane potential i.e. dendritic compartment 
du_bar_minus_dt(u_bar_minus, v) = (-u_bar_minus .+ v) ./ tau_u_bar_minus
du_bar_plus_dt(u_bar_plus, v) = (-u_bar_plus .+ v) ./ tau_u_bar_plus

## Normal STDP 
aplus = 0.2*nS; tplus = 20*ms; aminus = 0.25*nS; tminus = 20*ms
function stdp(t)
    if t > 0 
        return aplus * exp(abs(t)/tplus)
    else 
        return -aminus * exp(-abs(t)/tminus)
    end 
end 

A_LTD = 14e-5; A_LTP = 8e-5; theta_plus = -45.3*mV; theta_minus = -70.6*mV 
## Voltage-based STDP 
## X is the presynaptic spike train 
function vb_stdp(X, u, u_bar_plus, u_bar_minus, x_bar)
    new_weight = -A_LTD * X * (u_bar_minus .- theta_minus) + A_LTP * x_bar * (u - theta_plus) * (u_bar_plus - theta_minus)
    return new_weight
end 

## TODO - STD: E->PV
## TODO - STP: E->SST 
u_ESST = zeros(nr_pyc, nr_sst)
R_ESST = zeros(nr_pyc, nr_sst)
u_EPV  = zeros(nr_pyc, nr_pv)
R_EPV  = zeros(nr_pyc, nr_pv)

u_SSTPV  = zeros(nr_sst, nr_pv)
R_SSTPV  = zeros(nr_sst, nr_pv)
u_PVSST  = zeros(nr_pv, nr_sst)
R_PVSST  = zeros(nr_pv, nr_sst)

U_stp = 0.2; tau_rec_stp = 125*ms; tau_fac_stp = 500*ms
U_std = 0.25; tau_rec_std = 700*ms; tau_fac_std = 20*ms

tau_syn = 5*ms

## Simulation
for t = 2:steps
    v_d[:, t], w_d[:, t], v_s[:, t], w_s[:, t], I_dbg[:, t], I_sbg[:, t], t_pyc[:, t], t_soma[:] = simulatePyC(t, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], I_inj_d[:, t-1], I_inj_s[:, t-1], I_dbg[:, t-1], I_sbg[:, t-1], t_pyc[:, t-1], t_soma, W_SSTEd, W_PVEs, st_SSTEd, st_PVEs)

    ## Simulate for each type of interneuron 
    v_sst[:, t], I_sstbg[:, t], t_sst[:, t], tspike_sst[:] = simulateI(t, v_sst[:, t-1], I_sstbg[:, t-1], t_pyc[:, t-1], tspike_sst, W_ESST, W_PVSST, st_EsSST, st_PVSST) 
    v_pv[:, t], I_pvbg[:, t], t_pv[:, t], tspike_pv[:] = simulateI(t, v_pv[:, t-1], I_pvbg[:, t-1], t_pyc[:, t-1], tspike_pv, W_EPV, W_SSTPV, st_EsPV, st_SSTPV) 

    u_bar_minus[:] += du_bar_minus_dt(u_bar_minus, v_d[:, t]) .* dt
    u_bar_plus[:] += du_bar_plus_dt(u_bar_plus, v_d[:, t]) .* dt

    ## Update synaptic trace 
    st_EsSST[:] += (- st_EsSST ./ tau_syn + t_pyc[:, t]) .* dt
    st_EsPV[:] += (- st_EsPV ./ tau_syn + t_pyc[:, t]) .* dt
    st_SSTEd[:] += (- st_SSTEd ./ tau_syn + t_sst[:, t]) .* dt
    st_PVEs[:] += (- st_PVEs ./ tau_syn + t_pv[:, t]) .* dt
    st_SSTPV[:] += (- st_SSTPV ./ tau_syn + t_sst[:, t]) .* dt
    st_PVSST[:] += (- st_PVSST ./ tau_syn + t_pv[:, t]) .* dt

    ## Update weights - according various plasticity rules 
    ## E->I{SST,PV}
    ## E->SST - STP 
    for i=1:nr_pyc, j=1:nr_sst
        u_ESST[i, j] += ((U_stp - u_ESST[i, j]) / tau_fac_stp + U_stp * (1 - u_ESST[i, j]) * t_pyc[i, t]) * dt
        R_ESST[i, j] += ((1 - R_ESST[i, j]) / tau_rec_stp - u_ESST[i, j] * R_ESST[i, j] * t_pyc[i, t]) * dt
    
        W_ESST[i, j] = u_ESST[i, j] * R_ESST[i, j]
    end 
    ## E->PV - STD
    for i=1:nr_pyc, j=1:nr_pv
        u_EPV[i, j] += ((U_std - u_EPV[i, j]) / tau_fac_std + U_std * (1 - u_EPV[i, j]) * t_pyc[i, t]) * dt
        R_EPV[i, j] += ((1 - R_EPV[i, j]) / tau_rec_std - u_EPV[i, j] * R_EPV[i, j] * t_pyc[i, t]) * dt
    
        W_EPV[i, j] = u_EPV[i, j] * R_EPV[i, j]
    end 

    ## I{SST,PV}->E{s,d}
    for i=1:nr_sst, j=1:nr_pyc
        W_SSTEd[i, j] += vb_stdp(t_sst[i, t], v_d[j, t], u_bar_plus[j], u_bar_minus[j], st_SSTEd[i,j]) * dt
    end 
    for i=1:nr_pv, j=1:nr_pyc
        W_PVEs[i, j] += stdp(t_soma[j] - tspike_pv[i])
    end 

    ## I{SST,PV}->I{SST,PV}
    for i=1:nr_sst, j=1:nr_pv
        u_SSTPV[i, j] += ((U_std - u_SSTPV[i, j]) / tau_fac_std + U_std * (1 - u_SSTPV[i, j]) * t_sst[i, t]) * dt
        R_SSTPV[i, j] += ((1 - R_SSTPV[i, j]) / tau_rec_std - u_SSTPV[i, j] * R_SSTPV[i, j] * t_sst[i, t]) * dt
    
        W_SSTPV[i, j] = u_SSTPV[i, j] * R_SSTPV[i, j]
    end 
    for i=1:nr_pv, j=1:nr_sst
        u_PVSST[i, j] += ((U_std - u_PVSST[i, j]) / tau_fac_std + U_std * (1 - u_PVSST[i, j]) * t_pv[i, t]) * dt
        R_PVSST[i, j] += ((1 - R_PVSST[i, j]) / tau_rec_std - u_PVSST[i, j] * R_PVSST[i, j] * t_pv[i, t]) * dt

        W_PVSST[i, j] = u_PVSST[i, j] * R_PVSST[i, j]
    end 
end 

step_list = [1:steps;]
p1 = plot(step_list, v_s[:], label="Somatic Voltage")
p2 = plot(step_list, v_d[:], label="Dendritic Voltage")
p3 = plot(step_list, v_sst[:], label="SST Voltage")
p4 = plot(step_list, v_pv[:], label="PV Voltage")

p5 = plot(step_list, I_sbg[:], label="Somatic Background Current")
p6 = plot(step_list, I_dbg[:], label="Dendritic Background Current")
p7 = plot(step_list, I_sstbg[:], label="SST Background Current")
p8 = plot(step_list, I_pvbg[:], label="PV Background Current")

p9 = plot(step_list, t_pyc[:], label="PyC Spike Train")
p10 = plot(step_list, t_sst[:], label="SST Spike Train")
p11 = plot(step_list, t_pv[:], label="PV Spike Train") 

# Plot voltage trace
display(plot(p1, p3, p4, layout=(3,1)))

# Plot voltage trace and background current 
# display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(4,2), size=(1000,500)))

# Plot spike train 
# display(plot(p1, p2, p3, p4, p9, p10, p11))

end 