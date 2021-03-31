include("Simulation/NeuralNetwork.jl")

using Plots
using Formatting
using CSV, DataFrames
using .NeuralNetwork: Model.Connectivity.initialise_weights
using .NeuralNetwork: start_simulation
using .NeuralNetwork: Model.InjectedCurrent.initialise_inj_current, Model.InjectedCurrent.add_inj_current!

## Network size
nr_pyc = 100; nr_sst = 10; nr_pv = 10

## Simulation timing
t = 3000e-3; dt = 1e-3

## Stimulation
I_inj_amount = 10.2e-6 #nA
I_inj_duration = 30 #ms

## Selection

I_inj_d = zeros(nr_pyc, Int(t/dt)); I_inj_s = zeros(nr_pyc, Int(t/dt))

# add_inj_current!(I_inj_s, I_inj_amount, 1, 30, 24)
add_inj_current!(I_inj_s, I_inj_amount, 1, 30, 50)
# add_inj_current!(I_inj_s, I_inj_amount, 1, 30, 26)

add_inj_current!(I_inj_d, 10*I_inj_amount, 20, 30, 49)
add_inj_current!(I_inj_d, 10*I_inj_amount, 20, 30, 50)
add_inj_current!(I_inj_d, 10*I_inj_amount, 20, 30, 51)

# for kappa_excexc=14:2:20

# Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "alpha_excexc" => 2.6e-7,
                   "alpha_excinh" => 1.9e-7,
                   "alpha_inhexc" => 1.1e-7,
                   "alpha_inhinh" => 1e-7,
                   # Concentration parameter kappa
                   "kappa_excexc" => 18,
                   "kappa_excinh" => 0.01,
                   "kappa_inhexc" => 0.01,
                   "kappa_inhinh" => 0.01)

W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(nr_pyc, nr_sst, nr_pv, con_params)
v_s = start_simulation(t, dt, nr_pyc, nr_pv, nr_sst, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

# Write to CSV file
# filename = format("{1:d}.csv", kappa_excexc)
# CSV.write("C:/Users/Orion/Documents/University/Dissertation/Julia/data/volt-kappa_excexc/$filename", DataFrame(v_s), header=false)

# Plot somatic voltage
Plots.heatmap(v_s, c=:balance)

# end
