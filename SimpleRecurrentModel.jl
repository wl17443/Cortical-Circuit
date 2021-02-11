## Simple Model with PyC and INs 
module SimpleRecurrentNetwork 

include("PyramidalNeuron.jl")
include("Interneuron.jl")
include("Connectivity.jl")
include("Units.jl")
include("ModellingParameters.jl")

using .PyramidalNeuron
using .Interneuron
using .Connectivity
using .Units    
using .ModellingParameters
using Random, Distributions

## Variables 
## Voltages 
v_d = zeros(nr_pyc, steps)
w_d = zeros(nr_pyc, steps)
v_s = zeros(nr_pyc, steps)
w_s = zeros(nr_pyc, steps)
v_sst = zeros(nr_sst, steps)
v_pv = zeros(nr_pv, steps)
v_chc = zeros(nr_chc, steps)

I_sbg = randn(nr_pyc)
I_dbg = randn(nr_pyc)
I_sstbg = randn(nr_sst)
I_pvbg = randn(nr_pv)

## Initial values 
v_d[1] = -70*mV
v_s[1] = -70*mV
w_d[:, 1] = randn(nr_pyc)
w_s[:, 1] = randn(nr_pyc)
v_sst[1] = -70*mV
v_pv[1] = -70*mV

## Utilisation variable and recovery variable for Interneurons
u_sst = randn(nr_sst)
u_pv = randn(nr_pv)
R_sst = randn(nr_sst)
R_pv = randn(nr_pv)

## Spike timing
t_pyc = [[-1.0*s] for i=1:nr_pyc]
t_sst = [[-1.0*s] for i=1:nr_sst]
t_pv = [[-1.0*s] for i=1:nr_pv]

## Synaptic trace 
st_EsSST = randn(nr_sst)
st_EsPV = randn(nr_pv)
st_SSTEd = randn(nr_pyc)
st_PVEs = randn(nr_pyc)

## Simulation
for t = 2:steps
    v_dPrime, w_dPrime, v_sPrime, w_sPrime, I_dbgPrime, I_sbgPrime, t_Prime = simulatePyC(t, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], I_dbg, I_sbg, t_pyc, st_SSTEd, st_PVEs)

    ## Simulate for each type of interneuron 
    v_sstPrime, I_sstbgPrime = simulateI(t, v_sst[:, t-1], I_sstbg, u_sst, R_sst, W_ESST, W_PVSST) 
    v_pvPrime, I_pvbgPrime = simulateI(t, v_pv[:, t-1], I_pvbg, u_pv, R_pv, W_EPV, W_SSTPV) 

    t_ = t_Prime 
    v_d[:, t] = v_dPrime
    w_d[:, t] = w_dPrime
    v_s[:, t] = v_sPrime
    w_s[:, t] = w_sPrime 

    v_sst[:, t] = v_sstPrime 
    v_pv[:, t] = v_pvPrime 

    # I_dbg = I_dbgPrime 
    # I_sbg = I_sbgPrime 
    # I_sstbg = I_sstbgPrime 
    # I_pvbg = I_pvbgPrime 
end 

end 
