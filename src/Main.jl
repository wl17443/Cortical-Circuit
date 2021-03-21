include("Simulation/NeuralNetwork.jl")

using CSV, DataFrames
using .NeuralNetwork: Model.Connectivity.initialise_weights
using .NeuralNetwork: start_simulation

## Network size
nr_pyc = 50; nr_sst = 50; nr_pv = 50

## Simulation timing
t = 2000e-3; dt = 1e-3

## Background noise level
bgnoise_lvl = 1

## Connectivity
alpha_exc = 1e-9:1e-9:20e-9
alpha_inh = 1e-9:1e-9:20e-9

## Injected Current
I_inj_amount = 5e-6
I_inj_duration = 100 #ms

I_inj_d = zeros(nr_pyc, Int(t/dt)); I_inj_s = zeros(nr_pyc, Int(t/dt))

I_inj_s[25, 200:(200 + I_inj_duration)] .= I_inj_amount

for exc in alpha_exc, inh in alpha_inh
    W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(nr_pyc, nr_sst, nr_pv, exc, inh)
    v_s = start_simulation(t, dt, nr_pyc, nr_pv, nr_sst, bgnoise_lvl, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)
    CSV.write("C:/Users/Orion/Documents/University/Dissertation/Julia/csv/$exc-$inh-$I_inj_amount-$I_inj_duration.csv", DataFrame(v_s), writeheader=false)
end
