using CSV
using DelimitedFiles

function start_simulation(log, t::Float64, dt::Float64, analysis_slice::Int, network_params::Dict, noise_lvl, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

steps = Int(t/dt)
tau_syn = 5e-3

## Variables
## Voltages
v_d = zeros(network_params["nr_pyc"], steps)
v_s = zeros(network_params["nr_pyc"], steps)
w_d = zeros(network_params["nr_pyc"], steps)
w_s = zeros(network_params["nr_pyc"], steps)
v_sst = zeros(network_params["nr_sst"], steps)
v_pv = zeros(network_params["nr_pv"], steps)

## External background currents
I_sbg = zeros(network_params["nr_pyc"], steps)
I_dbg = zeros(network_params["nr_pyc"], steps)
I_sstbg = zeros(network_params["nr_sst"], steps)
I_pvbg = zeros(network_params["nr_pv"], steps)

## Spike train
t_pyc = zeros(network_params["nr_pyc"], steps)
t_sst = zeros(network_params["nr_sst"], steps)
t_pv = zeros(network_params["nr_pv"], steps)

## Last spike timing
t_soma = zeros(network_params["nr_pyc"])
tspike_sst = zeros(network_params["nr_sst"])
tspike_pv = zeros(network_params["nr_pv"])

## Synaptic trace
st_ESST = zeros(network_params["nr_pyc"], network_params["nr_sst"],)
st_EPV = zeros(network_params["nr_pyc"], network_params["nr_pv"])
st_SSTE = zeros(network_params["nr_sst"], network_params["nr_pyc"]);
st_PVE = zeros(network_params["nr_pv"], network_params["nr_pyc"])
st_SSTPV = zeros(network_params["nr_sst"], network_params["nr_pv"])
st_PVSST = zeros(network_params["nr_pv"], network_params["nr_sst"])
st_EE = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

## Initial values
## Voltages
v_d[:, 1] .= -70e-3; v_s[:, 1] .= -70e-3
v_sst[:, 1] .= -70e-3; v_pv[:, 1] .= -70e-3

## Last analysis time point
last_checkpoint = 1
checkpoint = 1

## EI current distribution
distr = zeros(7, Int(steps/analysis_slice))
firing_rate = zeros(3, Int(steps/analysis_slice))
location = zeros(Int(steps/analysis_slice))

## Simulation
for t = 2:steps
    v_d[:, t], w_d[:, t], v_s[:, t], w_s[:, t], I_dbg[:, t], I_sbg[:, t], t_pyc[:, t], t_soma[:] = simulatePyC(t, dt, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], noise_lvl, I_inj_d[:, t-1], I_inj_s[:, t-1], I_dbg[:, t-1], I_sbg[:, t-1], t_pyc[:, t-1], t_soma, st_SSTE, st_PVE, st_EE, W_SSTE, W_PVE, W_EE)

    ## Simulate for each type of interneuron
    v_sst[:, t], I_sstbg[:, t], t_sst[:, t], tspike_sst[:] = simulateI(t, dt, v_sst[:, t-1], noise_lvl, I_sstbg[:, t-1], tspike_sst, st_ESST, st_PVSST, W_ESST, W_PVSST)
    v_pv[:, t], I_pvbg[:, t], t_pv[:, t], tspike_pv[:] = simulateI(t, dt, v_pv[:, t-1], noise_lvl, I_pvbg[:, t-1], tspike_pv, st_EPV, st_SSTPV, W_EPV, W_SSTPV)

    ## Update synaptic trace
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_sst"]
        st_ESST[i, j] += (- st_ESST[i, j] ./ tau_syn) * dt + t_pyc[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pv"]
        st_EPV[i, j] += (- st_EPV[i, j] ./ tau_syn) * dt + t_pyc[i, t]
    end
    for i=1:network_params["nr_sst"], j=1:network_params["nr_pyc"]
        st_SSTE[i, j] += (- st_SSTE[i, j] ./ tau_syn) * dt + t_sst[i, t]
    end
    for i=1:network_params["nr_pv"], j=1:network_params["nr_pyc"]
        st_PVE[i, j] += (- st_PVE[i, j] ./ tau_syn) * dt + t_pv[i, t]
    end
    for i=1:network_params["nr_sst"], j=1:network_params["nr_pv"]
        st_SSTPV[i, j] += (- st_SSTPV[i, j] ./ tau_syn) * dt + t_sst[i, t]
    end
    for i=1:network_params["nr_pv"], j=1:network_params["nr_sst"]
        st_PVSST[i, j] += (- st_PVSST[i, j] ./ tau_syn) * dt + t_pv[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pyc"]
        st_EE[i, j] += (- st_EE[i, j] ./ tau_syn) * dt + t_pyc[i, t]
    end

    ## Analyse current activity
    if t % analysis_slice == 0 && log != false
        spread, activity_site = bump_status(v_s[:, last_checkpoint:t], analysis_slice, network_params["nr_pyc"])

        # avg_fr_pyc = sum(t_pyc[:, last_checkpoint:t]) / (analysis_slice * 1e-3)
        # avg_fr_sst = sum(t_sst[:, last_checkpoint:t]) / (analysis_slice * 1e-3)
        # avg_fr_pv = sum(t_pv[:, last_checkpoint:t]) / (analysis_slice * 1e-3)
        #
        # firing_rate[:, checkpoint] = [avg_fr_pyc, avg_fr_pv, avg_fr_sst]
        #
        # distr[:, checkpoint] = [sum(W_EE)/50 * avg_fr_pyc, sum(W_ESST)/50 * avg_fr_pyc, sum(W_EPV)/50 * avg_fr_pyc, sum(W_SSTE)/5 * avg_fr_sst, sum(W_PVE)/5 * avg_fr_pv, sum(W_SSTPV)/5 * avg_fr_sst, sum(W_PVSST)/5 * avg_fr_pv]
        # distr[:, checkpoint] = distr[:, checkpoint] ./ sum(distr[:, checkpoint])
        # map!(x -> 100 * x, distr[:, checkpoint], distr[:, checkpoint])
        #
        # location[checkpoint] = spread > 30 ? activity_site : 0

        write(log, "Checkpoint $checkpoint:\n")
        write(log, "Bump at neuron $activity_site, spreading across $spread neurons.\n")

        last_checkpoint = t
        checkpoint += 1
    end
end

# error = sum(abs.(diff(location)))

# return v_s, distr, firing_rate, error
return v_s

end
