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
v_d = zeros(nr_pyc, steps)
w_d = zeros(nr_pyc, steps)
v_s = zeros(nr_pyc, steps)
w_s = zeros(nr_pyc, steps)
v_sst = zeros(nr_sst, steps)
v_pv = zeros(nr_pv, steps)
v_chc = zeros(nr_chc, steps)

I_sbg = zeros(nr_pyc, steps)
I_dbg = zeros(nr_pyc, steps)
I_sstbg = zeros(nr_sst, steps)
I_pvbg = zeros(nr_pv, steps)

## Initial values 

v_d[1] = -70*mV
v_s[1] = -70*mV
w_d[:, 1] = randn(nr_pyc)
w_s[:, 1] = randn(nr_pyc)
v_sst[1] = -70*mV
v_pv[1] = -70*mV

u_sst = randn(nr_sst)
u_pv = randn(nr_pv)
R_sst = randn(nr_sst)
R_pv = randn(nr_pv)

## Spike timing
t_pycd = [[] for i=1:nr_pyc]
t_pycs = [[] for i=1:nr_pyc]
t_sst = [[] for i=1:nr_sst]
t_pv = [[] for i=1:nr_pv]

## Synaptic trace 
st_EsSST = randn(nr_sst)
st_EsPV = randn(nr_pv)

## Simulation
for t_step = 2:steps
    simulatePyC(t_step, dt, v_d[:, t_step-1], w_d[:, t_step-1], v_s[:, t_step-1], w_s[:, t_step-1], I_dbg[:, t_step-1], I_sbg[:, t_step-1], syntrace_sst, syntrace_pv)

    ## Simulate for each type of interneuron 
    simulateI(t_step, dt, v_sst[:, t_step-1], I_sstbg, u_sst, R_sst) 
    simulateI(t_step, dt, v_pv[:, t_step-1], I_pvbg, u_pv, R_pv) 
end 

end 
