module TestSST

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

W_ESST = zeros(nr_pyc, nr_sst)
W_PVSST = zeros(nr_pv, nr_sst)
U_ESST = rand(Uniform(0.1, 0.25), (nr_pyc, nr_sst))

v_sst = zeros(nr_sst, steps)

I_sstbg = zeros(nr_sst, steps)

st_EsSST = zeros(nr_pyc, nr_sst)
st_PVSST = zeros(nr_pv, nr_sst)

t_sst = zeros(nr_sst, steps)
t_pyc = zeros(nr_pyc, steps)

st_EsSST = zeros(nr_pyc, nr_sst)
st_PVSST = zeros(nr_pv, nr_sst)

v_sst[1] = -70*mV

u_sst = randn(nr_pyc, nr_sst)
R_sst = randn(nr_pyc, nr_sst)

for t = 2:steps
    ## Simulate for each type of interneuron 
    v_sst[:, t], I_sstbg[:, t], u_sst[:], R_sst[:], t_sst[:, t] = simulateI(t, v_sst[:, t-1], I_sstbg[:, t-1], t_pyc[:, t-1], u_sst, R_sst, U_ESST, W_ESST, W_PVSST, st_EsSST, st_PVSST) 

    ## Update synaptic trace 
    # st_EsSST[:] += (- st_EsSST ./ t_syn + t_pyc[:, t]) .* dt
    # st_PVSST[:] += (- st_PVSST ./ t_syn + t_pv[:, t]) .* dt
end 

step_list = [1:steps;]
p1 = plot(step_list, v_sst[:], label="SST Voltage")
display(plot(p1))

end 
