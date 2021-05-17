using CSV
using DelimitedFiles

function start_simulation(log, t, dt, analysis_slice, network_params, noise_lvl, Isinj, Idinj, W_EE, W_ESST, W_EPV, W_EVip, W_SSTE, W_PVE, W_PVVip, W_VipSST)

steps = Int(t/dt)
tau_syn = 5e-3

## Variables
## Voltages
v_d = zeros(network_params["nr_pyc"], steps)
v_s = zeros(network_params["nr_pyc"], steps)
v_sst = zeros(network_params["nr_sst"], steps)
v_pv = zeros(network_params["nr_pv"], steps)
v_vip = zeros(network_params["nr_vip"], steps)

## Activation variables
w_d = zeros(network_params["nr_pyc"], steps)
w_s = zeros(network_params["nr_pyc"], steps)

## External background currents
I_sbg = zeros(network_params["nr_pyc"])
I_dbg = zeros(network_params["nr_pyc"])
I_sstbg = zeros(network_params["nr_sst"])
I_pvbg = zeros(network_params["nr_pv"])
I_vipbg = zeros(network_params["nr_vip"])

## Spike train
t_pyc = zeros(network_params["nr_pyc"], steps)
t_sst = zeros(network_params["nr_sst"], steps)
t_pv = zeros(network_params["nr_pv"], steps)
t_vip = zeros(network_params["nr_vip"], steps)

## Spike timing of soma
t_soma = zeros(network_params["nr_pyc"])

## Synaptic trace
st_ESST = zeros(network_params["nr_pyc"], network_params["nr_sst"])
st_EPV = zeros(network_params["nr_pyc"], network_params["nr_pv"])
st_EVip = zeros(network_params["nr_pyc"], network_params["nr_vip"])
st_EE = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

st_SSTE = zeros(network_params["nr_sst"], network_params["nr_pyc"])
st_PVE = zeros(network_params["nr_pv"], network_params["nr_pyc"])

st_PVVip = zeros(network_params["nr_pv"], network_params["nr_vip"])
st_VipSST = zeros(network_params["nr_vip"], network_params["nr_sst"])

## Recovery variables u (Interneurons)
u_pv = zeros(network_params["nr_pv"])
u_sst = zeros(network_params["nr_sst"])
u_vip = zeros(network_params["nr_vip"])

## Initial values
## Voltages
v_d[:, 1] .= -60e-3; v_s[:, 1] .= -60e-3
v_sst[:, 1] .= -56; v_pv[:, 1] .= -55
v_vip[:, 1] .= -75

## Last analysis time point
last_checkpoint = 1
checkpoint = 1

## Simulation
for t = 2:steps
    v_d[:, t], w_d[:, t], v_s[:, t], w_s[:, t], I_dbg[:], I_sbg[:], t_pyc[:, t], t_soma[:] = simulatePyC(t, dt*1e-3, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], noise_lvl, Idinj[:, t-1], Isinj[:, t-1], I_dbg, I_sbg, t_pyc[:, t-1], t_soma, st_SSTE, st_PVE, st_EE, W_SSTE, W_PVE, W_EE)

    ## Simulate for each type of interneuron
    v_sst[:, t], u_sst[:], I_sstbg[:], t_sst[:, t] = simulateSST(dt, v_sst[:, t-1], u_sst, noise_lvl, I_sstbg, st_ESST, st_VipSST, W_ESST, W_VipSST)
    v_pv[:, t], u_pv[:], I_pvbg[:], t_pv[:, t] = simulatePV(dt, v_pv[:, t-1], u_pv, noise_lvl, I_pvbg, st_EPV, W_EPV)
    v_vip[:, t], u_vip[:], I_vipbg[:], t_vip[:, t] = simulateVip(dt, v_vip[:, t-1], u_vip, noise_lvl, I_vipbg, st_EVip, st_PVVip, W_EVip, W_PVVip)

    ## Update synaptic trace
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_sst"]
        st_ESST[i, j] += (- st_ESST[i, j] / tau_syn) * dt*1e-3 + t_pyc[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pv"]
        st_EPV[i, j] += (- st_EPV[i, j] / tau_syn) * dt*1e-3 + t_pyc[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_vip"]
        st_EVip[i, j] += (- st_EVip[i, j] / tau_syn) * dt*1e-3 + t_pyc[i, t]
    end
    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pyc"]
        st_EE[i, j] += (- st_EE[i, j] / tau_syn) * dt*1e-3 + t_pyc[i, t]
    end

    for i=1:network_params["nr_sst"], j=1:network_params["nr_pyc"]
        st_SSTE[i, j] += (- st_SSTE[i, j] / tau_syn) * dt*1e-3 + t_sst[i, t]
    end
    for i=1:network_params["nr_pv"], j=1:network_params["nr_pyc"]
        st_PVE[i, j] += (- st_PVE[i, j] / tau_syn) * dt*1e-3 + t_pv[i, t]
    end

    for i=1:network_params["nr_pv"], j=1:network_params["nr_vip"]
        st_PVVip[i, j] += (- st_PVVip[i, j] / tau_syn) * dt*1e-3 + t_pv[i, t]
    end
    for i=1:network_params["nr_vip"], j=1:network_params["nr_sst"]
        st_VipSST[i, j] += (- st_VipSST[i, j] / tau_syn) * dt*1e-3 + t_vip[i, t]
    end

    ## Analyse current activity
    if t % analysis_slice == 0 && log != false
        spread, activity_site = bump_status(v_s[:, last_checkpoint:t], analysis_slice, network_params["nr_pyc"])

        write(log, "Checkpoint $checkpoint:\n")
        write(log, "Bump at neuron $activity_site, spreading across $spread neurons.\n")

        last_checkpoint = t
        checkpoint += 1
    end
end

return v_s, v_d, v_pv, v_sst, v_vip

end
