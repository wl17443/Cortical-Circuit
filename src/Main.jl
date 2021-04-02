include("Simulation/NeuralNetwork.jl")

using Plots
using Dates
using Formatting
using CSV, DataFrames
using .NeuralNetwork: Model.Connectivity.initialise_weights
using .NeuralNetwork: start_simulation
using .NeuralNetwork: Model.InjectedCurrent.initialise_inj_current, Model.InjectedCurrent.add_inj_current!

## Network size
nr_pyc = 50; nr_sst = 5; nr_pv = 5

## Simulation timing
t = 6000e-3; dt = 1e-3

## Stimulation
I_inj_amount = 5e-6 #nA
I_inj_duration = 50 #ms

## Selection

I_inj_d = zeros(nr_pyc, Int(t/dt)); I_inj_s = zeros(nr_pyc, Int(t/dt))

add_inj_current!(I_inj_s, I_inj_amount, 1, I_inj_duration, 20)
# add_inj_current!(I_inj_s, I_inj_amount, 1, I_inj_duration, 40)

add_inj_current!(I_inj_d, 10*I_inj_amount, 25, I_inj_duration, 22)

# for alpha_excexc=3e-7:1e-8:4e-7, alpha_excinh=1e-7:1e-7:3e-7, alpha_inhexc=5e-7:2e-7:27e-7, kappa_excexc=16:2:30

# Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "alpha_excexc" => 6e-7,
                   "alpha_excinh" => 45e-7,
                   "alpha_inhexc" => 200e-7,
                   "alpha_inhinh" => 0.1e-7,
                   # Concentration parameter kappa
                   "kappa_excexc" => 16)

W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(nr_pyc, nr_sst, nr_pv, con_params)
v_s = start_simulation(t, dt, nr_pyc, nr_pv, nr_sst, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

# Write to CSV file
datetime = Dates.format(now(), "dd-mm-yy-HH-MM-SS")
filename = datetime * ".csv"
CSV.write("C:/Users/Orion/Documents/University/Dissertation/Julia/data/01-04-21-13-41/$filename", DataFrame(v_s), header=false)

# Plot somatic voltage
Plots.heatmap(v_s, c=:balance)

# end
