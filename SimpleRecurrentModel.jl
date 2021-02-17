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
using Noise 

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
 
I_sbg = zeros(nr_pyc, steps)
I_dbg = zeros(nr_pyc, steps)
I_sstbg = zeros(nr_sst, steps)
I_pvbg = zeros(nr_pv, steps)

## Initial values 
v_d[1] = -70*mV
v_s[1] = -70*mV
w_d[:, 1] = randn(nr_pyc)*nA
w_s[:, 1] = randn(nr_pyc)*nA

v_sst[1] = -70*mV
v_pv[1] = -70*mV

## Utilisation variable and recovery variable for Interneurons
u_sst = randn(nr_sst)
u_pv = randn(nr_pv)
R_sst = randn(nr_sst)
R_pv = randn(nr_pv)

## Last spike time
t_pyc = zeros(nr_pyc) 
t_sst = zeros(nr_sst)
t_pv = zeros(nr_pv)

## Synaptic trace 
st_EsSST = zeros(nr_sst)
st_EsPV = zeros(nr_pv)
st_SSTEd = zeros(nr_pyc)
st_PVEs = zeros(nr_pyc)
st_SSTPV = zeros(nr_pv)
st_PVSST = zeros(nr_sst)

I_dbg[:, 1] = add_gauss(-300*pA/2*ms .* ones(nr_pyc), 450*pA) #map!( x -> x < 0 ? 0 : x, I_dbg[:, 1], I_dbg[:, 1])
I_sbg[:, 1] = add_gauss(400*pA/2*ms .* ones(nr_pyc), 450*pA) #map!( x -> x < 0 ? 0 : x, I_sbg[:, 1], I_sbg[:, 1])
I_sstbg[:, 1] = add_gauss(-100*pA/2*ms .* ones(nr_pyc), 400*pA) #map!( x -> x < 0 ? 0 : x, I_sstbg[:, 1], I_sstbg[:, 1])
I_pvbg[:, 1] = add_gauss(-100*pA/2*ms .* ones(nr_pyc), 400*pA) #map!( x -> x < 0 ? 0 : x, I_pvbg[:, 1], I_pvbg[:, 1])

I_inj_d = ones(nr_pyc, steps) .* 2000*pA
I_inj_s = ones(nr_pyc, steps) .* 2000*pA

## Simulation
for t = 2:steps
    v_d[:, t], w_d[:, t], v_s[:, t], w_s[:, t], I_dbg[:, t], I_sbg[:, t], t_pyc[:], st_EsSST[:], st_EsPV[:] = simulatePyC(t, v_d[:, t-1], w_d[:, t-1], v_s[:, t-1], w_s[:, t-1], I_inj_d[:, t-1], I_inj_s[:, t-1], I_dbg[:, t-1], I_sbg[:, t-1], t_pyc, W_SSTEd, W_PVEs, st_SSTEd, st_PVEs, st_EsSST, st_EsPV)

    ## Simulate for each type of interneuron 
    v_sst[:, t], I_sstbg[:, t], u_sst[:], R_sst[:], t_sst[:], st_SSTEd[:], st_SSTPV[:] = simulateI(t, v_sst[:, t-1], I_sstbg[:, t-1], t_sst, u_sst, R_sst, U_ESST, W_ESST, W_PVSST, st_EsSST, st_PVSST, st_SSTEd, st_SSTPV) 
    v_pv[:, t], I_pvbg[:, t], u_pv[:], R_pv[:], t_pv[:], st_PVEs[:], st_PVSST[:] = simulateI(t, v_pv[:, t-1], I_pvbg[:, t-1], t_pv, u_pv, R_pv, U_EPV, W_EPV, W_SSTPV, st_EsPV, st_SSTPV, st_PVEs, st_PVSST) 
end 

step_list = [1:steps;]
p1 = plot(step_list, v_s[:], label="Somatic Voltage")
p2 = plot(step_list, v_d[:], label="Dendritic Voltage")
p3 = plot(step_list, v_sst[:], label="SST Voltage")
p4 = plot(step_list, v_pv[:], label="PV Voltage")
p5 = plot(step_list, I_sbg[:], label="Somatic bg Current")
p6 = plot(step_list, I_dbg[:], label="Dendritic bg Current")
p7 = plot(step_list, I_sstbg[:], label="SST bg Current")
p8 = plot(step_list, I_pvbg[:], label="PV bg Current")

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(4,2), size=(1000,1000)))

end 
