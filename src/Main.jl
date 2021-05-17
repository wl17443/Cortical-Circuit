include("Simulation/NeuralNetwork.jl")
include("Analysis/WriteData/WriteData.jl")
include("Analysis/Visualisation/Visualise.jl")

using CSV
using Plots
using Dates
using DelimitedFiles
using Random, Distributions
using .WriteData: write2csv
using .Visualise: save_heatmap, show_heatmap
using .NeuralNetwork: Model.Connectivity.initialise_weights
using .NeuralNetwork: start_simulation
using .NeuralNetwork: Model.InjectedCurrent.initialise_inj_current, Model.InjectedCurrent.add_inj_current!

network_params = Dict( # Network parameters - size
                       "nr_pyc" => 100,
                       "nr_sst" => 5,
                       "nr_pv"  => 5,
                       "nr_vip" => 10)

## Simulation timing
t = 5000; dt = 0.1 # ms

## Stimulation
I_inj_s = zeros(network_params["nr_pyc"], Int(t/dt))
stimulation_strength = 5e-6 # uA
stimulation_duration = 500 # 50 ms
add_inj_current!(I_inj_s, stimulation_strength, 1, stimulation_duration, 20)

## Selection
I_inj_d = zeros(network_params["nr_pyc"], Int(t/dt))
selection_strength = 10*stimulation_strength
selection_duration = 50 # 5 ms
add_inj_current!(I_inj_d, selection_strength, 10, selection_duration, 20)

noise = 0.0

# Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "s_excexc" => 3550e-9,
                   "s_excsst" => 0,
                   "s_excpv"  => 100,
                   "s_excvip" => 0,

                   "s_sstexc" => 0e-9,
                   "s_pvexc"  => 90e-9,

                   "s_pvvip"  => 0,
                   "s_vipsst"  => 0,

                   # Concentration parameter kappa
                   "k_excexc" => 80 )

noise_lvl = Dict( "som" => 1000,
                  "den" => 1000,
                  "inh" => 0)

W_EE, W_ESST, W_EPV, W_EVip, W_SSTE, W_PVE, W_PVVip, W_VipSST = initialise_weights(network_params, con_params, noise)

v_s, v_d, v_pv, v_sst, v_vip = start_simulation(false, t, dt, 500, network_params, noise_lvl, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_EVip, W_SSTE, W_PVE, W_PVVip, W_VipSST)

display(Plots.heatmap(v_s, title="Somatic Voltage Trace", xaxis="Time (ms)", yaxis="Nr. of Neuron", c=:balance))
# display(plot(collect(1:1000), v_pv[2, 1:1000]))
# Plots.heatmap(W_EE)
