include("Simulation/NeuralNetwork.jl")

using Plots
using Formatting
using CSV, DataFrames
using .NeuralNetwork: Model.Connectivity.initialise_weights
using .NeuralNetwork: start_simulation
using .NeuralNetwork: Model.InjectedCurrent.initialise_inj_current, Model.InjectedCurrent.add_inj_current!

## Network size
nr_pyc = 50; nr_sst = 50; nr_pv = 50

## Simulation timing
t = 6000e-3; dt = 1e-3

## Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "alpha_excexc" => 15e-9 ,
                   "alpha_excinh" => 22e-9,
                   "alpha_inhexc" => 1e-9,
                   "alpha_inhinh" => 0.94e-9,
                   # Concentration parameter kappa
                   "kappa_excexc" => 20,
                   "kappa_excinh" => 0.01,
                   "kappa_inhexc" => 0.01,
                   "kappa_inhinh" => 0.01)

## Injected Current
I_inj_amount = 6.2e-6 #nA
I_inj_duration = 50 #ms

I_inj_d = zeros(nr_pyc, Int(t/dt)); I_inj_s = zeros(nr_pyc, Int(t/dt))
add_inj_current!(I_inj_d, I_inj_amount, 1000, 50, 16)
add_inj_current!(I_inj_s, I_inj_amount, 1, I_inj_duration, 10)

W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(nr_pyc, nr_sst, nr_pv, con_params)
v_s, v_d, _, _ = start_simulation(t, dt, nr_pyc, nr_pv, nr_sst, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

# filename = format("{1:.1e}-{2:.1e}.csv", I_inj_amount, I_inj_duration)
# CSV.write("C:/Users/Orion/Documents/University/Dissertation/Julia/csv/two_activities/$filename", DataFrame(v_s), writeheader=false)

display(Plots.heatmap(v_s, c=:balance))
display(Plots.heatmap(v_d, c=:balance))
