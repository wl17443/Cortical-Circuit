## Main Network Model 
module RingAttractorNetwork

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

steps = Int(t/dt)

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
I_chcbg = zeros(nr_chc, steps)

## Initial values 
v_d[1] = -70*mV
w_d[1] = randn(nr_pyc)
v_s[1] = -70*mV
w_s[1] = randn(nr_pyc)

I_sbg[1] = 0*nA
I_dbg[1] = 0*nA

v_chc[1] = -70*mV
v_sst[1] = -70*mV
v_pv[1] = -70*mV

u_sst = randn(nr_sst)
u_pv = randn(nr_pv)
u_chc = randn(nr_chc)

R_sst = randn(nr_sst)
R_pv = randn(nr_pv)
R_chc = randn(nr_chc)

t_ = zeros(nr_pyc)

## Simulation 
for t_step = 2:(steps)
    simulatePyC(t_step, dt, v_d[:, t_step-1], w_d[:, t_step-1], v_s[:, t_step-1], w_s[:, t_step-1], I_dbg[:, t_step-1], I_sbg[:, t_step-1])
    
    ## Simulate for each type of interneuron 
    simulateI(t_step, dt, v_sst[:, t_step-1], I_sstbg, u_sst, R_sst) 
    simulateI(t_step, dt, v_pv[:, t_step-1], I_pvbg, u_pv, R_pv) 
    simulateI(t_step, dt, v_chc[:, t_step-1], I_chcbg, u_chc, R_chc) 
end 

# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 
