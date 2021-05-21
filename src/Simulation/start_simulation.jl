#=============================================================================================================

Run simulation with this function. 

=============================================================================================================#

using CSV
using DelimitedFiles

function start_simulation(t, dt, analysis_slice, network_params, Isinj, Idinj, W_EE, W_ESST, W_EPV, W_EVip, W_SSTE, W_PVE, W_PVVip, W_VipSST, D)

steps = Int(t/dt)
tau_syn = 5e-3

## Variables
## Voltages
v1_d = zeros(network_params["nr_pyc"], steps)
v1_s = zeros(network_params["nr_pyc"], steps)
v2_d = zeros(network_params["nr_pyc"], steps)
v2_s = zeros(network_params["nr_pyc"], steps)

v_sst = zeros(network_params["nr_sst"], steps)
v_pv = zeros(network_params["nr_pv"], steps)
v_vip = zeros(network_params["nr_vip"], steps)

## Activation variables
w1_d = zeros(network_params["nr_pyc"], steps)
w1_s = zeros(network_params["nr_pyc"], steps)
w2_d = zeros(network_params["nr_pyc"], steps)
w2_s = zeros(network_params["nr_pyc"], steps)

## Spike train
t1_pyc = zeros(network_params["nr_pyc"], steps)
t2_pyc = zeros(network_params["nr_pyc"], steps)

t_sst = zeros(network_params["nr_sst"], steps)
t_pv = zeros(network_params["nr_pv"], steps)
t_vip = zeros(network_params["nr_vip"], steps)

## Spike timing of soma
t1_soma = zeros(network_params["nr_pyc"])
t2_soma = zeros(network_params["nr_pyc"])

## Synaptic trace
st_E1SST = zeros(network_params["nr_pyc"], network_params["nr_sst"])
st_E1PV = zeros(network_params["nr_pyc"], network_params["nr_pv"])
st_E1Vip = zeros(network_params["nr_pyc"], network_params["nr_vip"])
st_E1E1 = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

st_E2SST = zeros(network_params["nr_pyc"], network_params["nr_sst"])
st_E2PV = zeros(network_params["nr_pyc"], network_params["nr_pv"])
st_E2Vip = zeros(network_params["nr_pyc"], network_params["nr_vip"])
st_E2E2 = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

st_E1E2 = zeros(network_params["nr_pyc"], network_params["nr_pyc"])
st_E2E1 = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

st_SSTE1 = zeros(network_params["nr_sst"], network_params["nr_pyc"])
st_PVE1 = zeros(network_params["nr_pv"], network_params["nr_pyc"])

st_SSTE2 = zeros(network_params["nr_sst"], network_params["nr_pyc"])
st_PVE2 = zeros(network_params["nr_pv"], network_params["nr_pyc"])

st_PVVip = zeros(network_params["nr_pv"], network_params["nr_vip"])
st_VipSST = zeros(network_params["nr_vip"], network_params["nr_sst"])

## Recovery variables u (Interneurons)
u_pv = zeros(network_params["nr_pv"])
u_sst = zeros(network_params["nr_sst"])
u_vip = zeros(network_params["nr_vip"])

## Firing rates
fr_pyc1= zeros(network_params["nr_pyc"], Int(steps/analysis_slice))
fr_pyc2= zeros(network_params["nr_pyc"], Int(steps/analysis_slice))
fr_pv= zeros(network_params["nr_pv"], Int(steps/analysis_slice))
fr_sst= zeros(network_params["nr_sst"], Int(steps/analysis_slice))
fr_vip= zeros(network_params["nr_vip"], Int(steps/analysis_slice))

## Initial values
## Voltages
v1_d[:, 1] .= -60e-3; v1_s[:, 1] .= -60e-3
v2_d[:, 1] .= -60e-3; v2_s[:, 1] .= -60e-3

v_sst[:, 1] .= -56; v_pv[:, 1] .= -55
v_vip[:, 1] .= -75

## Last analysis time point
last_checkpoint = 1
checkpoint = 1

mean_location1 = zeros(Int(steps/analysis_slice))
var_location1 = zeros(Int(steps/analysis_slice))
mean_location2 = zeros(Int(steps/analysis_slice))
var_location2 = zeros(Int(steps/analysis_slice))

## Simulation
for t = 2:steps
    v1_d[:, t], w1_d[:, t], v1_s[:, t], w1_s[:, t], t1_pyc[:, t], t1_soma[:] = simulatePyC1(t, dt*1e-3, v1_d[:, t-1], w1_d[:, t-1], v1_s[:, t-1], w1_s[:, t-1], Idinj[:, t-1], t1_pyc[:, t-1], t1_soma, st_SSTE1, st_PVE1, st_E1E1, st_E2E1, W_SSTE, W_PVE, W_EE, D)
    v2_d[:, t], w2_d[:, t], v2_s[:, t], w2_s[:, t], t2_pyc[:, t], t2_soma[:] = simulatePyC2(t, dt*1e-3, v2_d[:, t-1], w2_d[:, t-1], v2_s[:, t-1], w2_s[:, t-1], Isinj[:, t-1], t2_pyc[:, t-1], t2_soma, st_SSTE2, st_PVE2, st_E2E2, st_E1E2, W_SSTE, W_PVE, W_EE, D)

    ## Simulate for each type of interneuron
    v_sst[:, t], u_sst[:], t_sst[:, t] = simulateSST(dt, v_sst[:, t-1], u_sst, st_E1SST, st_E2SST, st_VipSST, W_ESST, W_VipSST, D)
    v_pv[:, t], u_pv[:], t_pv[:, t] = simulatePV(dt, v_pv[:, t-1], u_pv, st_E1PV, st_E2PV, W_EPV, D)
    v_vip[:, t], u_vip[:], t_vip[:, t] = simulateVip(dt, v_vip[:, t-1], u_vip, st_E1Vip, st_E2Vip, st_PVVip, W_EVip, W_PVVip, D)

    ## Update synaptic trace
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_sst"]
        st_E1SST[i, j] += (- st_E1SST[i, j] / tau_syn) * dt*1e-3 + t1_pyc[i, t]
        st_E2SST[i, j] += (- st_E2SST[i, j] / tau_syn) * dt*1e-3 + t2_pyc[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pv"]
        st_E1PV[i, j] += (- st_E1PV[i, j] / tau_syn) * dt*1e-3 + t1_pyc[i, t]
        st_E2PV[i, j] += (- st_E2PV[i, j] / tau_syn) * dt*1e-3 + t2_pyc[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_vip"]
        st_E1Vip[i, j] += (- st_E1Vip[i, j] / tau_syn) * dt*1e-3 + t1_pyc[i, t]
        st_E2Vip[i, j] += (- st_E2Vip[i, j] / tau_syn) * dt*1e-3 + t2_pyc[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pyc"]
        st_E1E1[i, j] += (- st_E1E1[i, j] / tau_syn) * dt*1e-3 + t1_pyc[i, t]
        st_E2E2[i, j] += (- st_E2E2[i, j] / tau_syn) * dt*1e-3 + t2_pyc[i, t]

        st_E1E2[i, j] += (- st_E1E2[i, j] / tau_syn) * dt*1e-3 + t1_pyc[i, t]
        st_E2E1[i, j] += (- st_E2E1[i, j] / tau_syn) * dt*1e-3 + t2_pyc[i, t]
    end

    for i=1:network_params["nr_sst"], j=1:network_params["nr_pyc"]
        st_SSTE1[i, j] += (- st_SSTE1[i, j] / tau_syn) * dt*1e-3 + t_sst[i, t]
        st_SSTE2[i, j] += (- st_SSTE2[i, j] / tau_syn) * dt*1e-3 + t_sst[i, t]
    end
    for i=1:network_params["nr_pv"], j=1:network_params["nr_pyc"]
        st_PVE1[i, j] += (- st_PVE1[i, j] / tau_syn) * dt*1e-3 + t_pv[i, t]
        st_PVE2[i, j] += (- st_PVE2[i, j] / tau_syn) * dt*1e-3 + t_pv[i, t]
    end

    for i=1:network_params["nr_pv"], j=1:network_params["nr_vip"]
        st_PVVip[i, j] += (- st_PVVip[i, j] / tau_syn) * dt*1e-3 + t_pv[i, t]
    end
    for i=1:network_params["nr_vip"], j=1:network_params["nr_sst"]
        st_VipSST[i, j] += (- st_VipSST[i, j] / tau_syn) * dt*1e-3 + t_vip[i, t]
    end

    ## Analyse current activity
    if t % analysis_slice == 0 # && log != false
        spread1, activity_site1 = bump_status(v1_s[:, last_checkpoint:t], analysis_slice, network_params["nr_pyc"])
        spread2, activity_site2 = bump_status(v2_s[:, last_checkpoint:t], analysis_slice, network_params["nr_pyc"])

        # fr_pyc1[:, checkpoint] = sum(t1_pyc[:, last_checkpoint:t], dims=2) ./ (analysis_slice*dt*1e-3)
        # fr_pyc2[:, checkpoint] = sum(t2_pyc[:, last_checkpoint:t], dims=2) ./ (analysis_slice*dt*1e-3)
        # fr_pv[:, checkpoint] = sum(t_pv[:, last_checkpoint:t], dims=2) ./ (analysis_slice*dt*1e-3)
        # fr_sst[:, checkpoint] = sum(t_sst[:, last_checkpoint:t], dims=2) ./ (analysis_slice*dt*1e-3)
        # fr_vip[:, checkpoint] = sum(t_vip[:, last_checkpoint:t], dims=2) ./ (analysis_slice*dt*1e-3)

        # write(log, "Checkpoint $checkpoint:\n")
        # write(log, "Bump at neuron $activity_site, spreading across $spread neurons.\n")

        mean_location1[checkpoint] = activity_site1
        var_location1[checkpoint] = spread1
        mean_location2[checkpoint] = activity_site2
        var_location2[checkpoint] = spread2

        last_checkpoint = t
        checkpoint += 1
    end
end

return v1_s, v2_s, v_pv, v_sst, v_vip
# return sum(abs.(diff(mean_location1))), count(x -> x >= 50, var_location1), sum(abs.(diff(mean_location2))), count(x -> x >= 50, var_location2)

end
