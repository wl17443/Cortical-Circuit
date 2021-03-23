include("Simulation/NeuralNetwork.jl")

using Plots
using Formatting
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
con_params = Dict( # Post synaptic current amplitude
                   "alpha_excexc" => 26e-9 ,
                   "alpha_excinh" => 21e-9,
                   "alpha_inhexc" => 2e-9,
                   "alpha_inhinh" => 1e-9,
                   # Concentration parameter kappa
                   "kappa_excexc" => 4,
                   "kappa_excinh" => 2,
                   "kappa_inhexc" => 0.5,
                   "kappa_inhinh" => 1)

## Injected Current
I_inj_amount = 5e-6 #nA
I_inj_duration = 50 #ms

I_inj_d = zeros(nr_pyc, Int(t/dt)); I_inj_s = zeros(nr_pyc, Int(t/dt))
I_inj_s[25, 200:(200+I_inj_duration)] .= I_inj_amount

W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(nr_pyc, nr_sst, nr_pv, con_params)
v_s = start_simulation(t, dt, nr_pyc, nr_pv, nr_sst, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

# filename = format("{1:.1e}-{2:.1e}-{3:.1e}-{4:d}.csv", exc, inh, I_inj_amount, I_inj_duration)
# CSV.write("C:/Users/Orion/Documents/University/Dissertation/Julia/csv/$filename", DataFrame(v_s), writeheader=false)

display(Plots.heatmap(v_s, c=:balance))
