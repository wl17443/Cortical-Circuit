include("../src/model_components/Units.jl")
include("../src/model_components/Connectivity.jl")
include("../src/network_models/RingAttractorModel.jl")

using .Units
using .Connectivity: initialise_weights
using .RingAttractorModel: start_simulation

nr_pyc = 50; nr_sst = 50; nr_pv = 50

## Simulation
t = 2000*ms; dt = 1*ms

## Background noise level
bgnoise_lvl = 0

## Connectivity
alpha_exc = 1*pA
alpha_inh = 5*pA

## Injected Current
I_inj_d = zeros(nr_pyc, Int(t/dt)); I_inj_s = zeros(nr_pyc, Int(t/dt))

I_inj_s[24, 200:300] .= 5*nA
I_inj_s[25, 200:300] .= 5*nA
I_inj_s[26, 200:300] .= 5*nA

W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(nr_pyc, nr_sst, nr_pv, alpha_exc, alpha_inh)

start_simulation(t, dt, nr_pyc, nr_pv, nr_sst, bgnoise_lvl, alpha_exc, alpha_inh, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)
