module TestPV

include("Units.jl")
include("Interneuron.jl")

using .Units 
using .Interneuron
using Random, Distributions 
using Plots 

nr_sst = 1 
nr_pyc = 1
nr_pv = 1 

t = 1000*ms 
dt = 1*ms 
steps = Int(t/dt)

W_EPV = zeros(nr_pyc, nr_pv)
W_SSTPV = zeros(nr_sst, nr_pv)
U_EPV = rand(Uniform(0.1, 0.25), (nr_pyc, nr_pv))

v_pv = zeros(nr_pv, steps)

I_pvbg = zeros(nr_pv, steps)

st_EsPV = zeros(nr_pyc, nr_pv)
st_SSTPV = zeros(nr_sst, nr_pv)

t_pv = zeros(nr_pv, steps)
t_pyc = zeros(nr_pyc, steps)

v_pv[1] = -70*mV

u_pv = randn(nr_pyc, nr_pv)
R_pv = randn(nr_pyc, nr_pv)

for t = 2:steps
    ## Simulate for each type of interneuron 
    v_pv[:, t], I_pvbg[:, t], u_pv[:], R_pv[:], t_pv[:, t] = simulateI(t, v_pv[:, t-1], I_pvbg[:, t-1], t_pyc[:, t-1], u_pv, R_pv, U_EPV, W_EPV, W_SSTPV, st_EsPV, st_SSTPV) 

    ## Update synaptic trace 
    # st_EsSST[:] += (- st_EsSST ./ t_syn + t_pyc[:, t]) .* dt
    # st_PVSST[:] += (- st_PVSST ./ t_syn + t_pv[:, t]) .* dt
end 

step_list = [1:steps;]
p1 = plot(step_list, v_pv[:], label="PV Voltage")
display(plot(p1))

end 
