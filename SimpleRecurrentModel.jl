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

## Weights being optimised 
W_ESST = randn((nr_pyc, nr_sst))/(nr_pyc*nr_sst)
W_EPV = randn((nr_pyc, nr_pv))/(nr_pyc*nr_pv)

W_SSTEd = randn((nr_sst, nr_pyc))*0.2/(nr_sst*nr_pyc)
W_PVEs = randn((nr_pv, nr_pyc))*0.2/(nr_pv*nr_pyc)

## Inter-Interneuron connectivity 
W_SSTPV = randn((nr_sst, nr_pv))/(nr_sst*nr_pv)
W_PVSST = randn((nr_pv, nr_sst))/(nr_pv*nr_sst)

## U - Initial release probability 
U_ESST = rand(Uniform(0.1, 0.25), (nr_pyc, nr_sst))
U_EPV = rand(Uniform(0.1, 0.25), (nr_pyc, nr_pv))

## Variables 
## Voltages 
v_d = zeros(nr_pyc, steps)
w_d = zeros(nr_pyc, steps)
v_s = zeros(nr_pyc, steps)
w_s = zeros(nr_pyc, steps)
v_sst = zeros(nr_sst, steps)
v_pv = zeros(nr_pv, steps)
# v_chc = zeros(nr_chc, steps)

I_sbg = ones(nr_pyc, steps)*nA
I_dbg = ones(nr_pyc, steps)*nA
I_sstbg = ones(nr_sst, steps)*nA
I_pvbg = ones(nr_pv, steps)*nA

## Initial values 
v_d[1] = -70*mV
v_s[1] = -70*mV
w_d[:, 1] = randn(nr_pyc)*nA
w_s[:, 1] = randn(nr_pyc)*nA

v_sst[1] = -70*mV
v_pv[1] = -70*mV

## Utilisation variable and recovery variable for Interneurons
u_sst = zeros(nr_sst)
u_pv = zeros(nr_pv)
R_sst = zeros(nr_sst)
R_pv = zeros(nr_pv)

## Last spike time
t_pyc = ones(nr_pyc) .* -1
t_sst = ones(nr_sst) .* -1
t_pv = ones(nr_pv) .* -1

## Synaptic trace 
st_EsSST = zeros(nr_sst)
st_EsPV = zeros(nr_pv)
st_SSTEd = zeros(nr_pyc)
st_PVEs = zeros(nr_pyc)
st_SSTPV = zeros(nr_pv)
st_PVSST = zeros(nr_sst)

I_inj_d = ones(nr_pyc, steps) * 1000*pA
I_inj_s = ones(nr_pyc, steps) * 1000*pA

## Simulation
for t = 2:steps
    v_dPrime, w_dPrime, v_sPrime, w_sPrime, I_dbgPrime, I_sbgPrime, t_pycPrime, st_EsSSTPrime, st_EsPVPrime = simulatePyC(t, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], I_inj_d[:, t-1], I_inj_s[:, t-1], I_dbg[:, t-1], I_sbg[:, t-1], t_pyc, W_SSTEd, W_PVEs, st_SSTEd, st_PVEs, st_EsSST, st_EsPV)

    ## Simulate for each type of interneuron 
    v_sstPrime, I_sstbgPrime, t_sstPrime, st_SSTEdPrime, st_SSTPVPrime = simulateI(t, v_sst[:, t-1], I_sstbg[:, t-1], t_sst, u_sst, R_sst, U_ESST, W_ESST, W_PVSST, st_EsSST, st_PVSST, st_SSTEd, st_SSTPV) 
    v_pvPrime, I_pvbgPrime, t_pvPrime, st_PVEsPrime, st_PVSSTPrime = simulateI(t, v_pv[:, t-1], I_pvbg[:, t-1], t_pv, u_pv, R_pv, U_EPV, W_EPV, W_SSTPV, st_EsPV, st_SSTPV, st_PVEs, st_PVSST) 

    t_pyc[:] = t_pycPrime 
    t_sst[:] = t_sstPrime 
    t_pv[:] = t_pvPrime 

    v_d[:, t] = v_dPrime
    w_d[:, t] = w_dPrime
    v_s[:, t] = v_sPrime
    w_s[:, t] = w_sPrime 

    v_sst[:, t] = v_sstPrime 
    v_pv[:, t] = v_pvPrime 

    I_dbg[:, t] = I_dbgPrime 
    I_sbg[:, t] = I_sbgPrime 
    I_sstbg[:, t] = I_sstbgPrime 
    I_pvbg[:, t] = I_pvbgPrime 

    st_EsPV[:] = st_EsPVPrime 
    st_EsSST[:] = st_EsSSTPrime
    st_SSTEd[:] = st_SSTEdPrime
    st_SSTPV[:] = st_SSTPVPrime
    st_PVEs[:] = st_PVEsPrime
    st_PVSST[:] = st_PVSSTPrime 
end 

step_list = [1:steps;]
p1 = plot(step_list, v_s[:], label="Somatic Voltage")
p2 = plot(step_list, v_d[:], label="Dendritic Voltage")
p3 = plot(step_list, v_sst[:], label="SST Voltage")
p4 = plot(step_list, v_pv[:], label="PV Voltage")
p5 = plot(step_list, I_dbg[:], label="Dedritic bg Current")
p6 = plot(step_list, I_sbg[:], label="Somatic bg Current")
p7 = plot(step_list, I_sstbg[:], label="SST bg Current")
p8 = plot(step_list, I_pvbg[:], label="PV bg Current")

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(4,2)))
end 
