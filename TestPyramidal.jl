module TestPyramidal 

include("Units.jl")
include("PyramidalNeuron.jl")

using .Units
using .PyramidalNeuron
using Random, Distributions 
using Plots 

nr_sst = 1 
nr_pyc = 1
nr_pv = 1 

t = 1000*ms 
dt = 1*ms 
steps = Int(t/dt)

W_SSTEd = zeros(nr_sst, nr_pyc)
W_PVEs = zeros(nr_pv, nr_pyc)

v_d = zeros(nr_pyc, steps)
v_s = zeros(nr_pyc, steps)
w_d = zeros(nr_pyc, steps)
w_s = zeros(nr_pyc, steps)

I_sbg = zeros(nr_pyc, steps)
I_dbg = zeros(nr_pyc, steps)

t_pyc = zeros(nr_pyc, steps) 

st_SSTEd = zeros(nr_sst, nr_pyc)
st_PVEs = zeros(nr_pv, nr_pyc)

v_d[1] = -70*mV
v_s[1] = -70*mV
w_d[:, 1] = randn(nr_pyc)*nA
w_s[:, 1] = randn(nr_pyc)*nA

I_inj_d = ones(nr_pyc, steps) .* 0*pA
I_inj_s = ones(nr_pyc, steps) .* 0*pA

## Simulation
for t = 2:steps
    v_d[:, t], w_d[:, t], v_s[:, t], w_s[:, t], I_dbg[:, t], I_sbg[:, t], t_pyc[:, t] = simulatePyC(t, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], I_inj_d[:, t-1], I_inj_s[:, t-1], I_dbg[:, t-1], I_sbg[:, t-1], t_pyc[:, t-1], W_SSTEd, W_PVEs, st_SSTEd, st_PVEs)
end 

step_list = [1:steps;]
p1 = plot(step_list, v_s[:], label="Somatic Voltage")
p2 = plot(step_list, v_d[:], label="Dendritic Voltage")

display(plot(p1, p2, layout=(2,1)))

end 
