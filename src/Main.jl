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
                       "nr_sst" => 10,
                       "nr_pv"  => 10,
                       "nr_vip" => 20 )

## Simulation timing
t = 3000 # ms
dt = 0.1 # ms

## Initialise injected current
## Stimulation
I_inj_s = zeros(network_params["nr_pyc"], Int(t/dt))
stimulation_strength = 5e-6 # uA
stimulation_duration = 500 # 50 ms
add_inj_current!(I_inj_s, stimulation_strength, 1, stimulation_duration, 20)
# add_inj_current!(I_inj_s, stimulation_strength, 1, stimulation_duration, 54)

## Selection
I_inj_d = zeros(network_params["nr_pyc"], Int(t/dt))
selection_strength = 10*stimulation_strength
selection_duration = 50 # 5 ms
# add_inj_current!(I_inj_d, selection_strength, 500, selection_duration, 20)
# add_inj_current!(I_inj_d, selection_strength, 4010, selection_duration, 40)

## Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "s_excexc"  => 3000e-9,
                   "s_excsst"  => 10, # 80,
                   "s_excpv"   => 40, # 20,
                   "s_excvip"  => 90,

                   "s_sstexc"  => 1000e-9, # 400e-9, # 100e-9,
                   "s_pvexc"   => 100e-9, # 50e-9,

                   "s_pvvip"   => 20,
                   "s_vipsst"  => 200, # 20,

                   # Concentration parameter kappa
                   "k_excexc" => 280 )

## Variance in same-type synaptic weights
percentage = 0.0

## Initialise synaptic weights
W_EE, W_ESST, W_EPV, W_EVip, W_SSTE, W_PVE, W_PVVip, W_VipSST = initialise_weights(network_params, con_params, percentage)

## Background noise variable
D = 1 * 450e-9
v1, v2, v_pv, v_sst, v_vip= start_simulation(t, dt, 500, network_params, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_EVip, W_SSTE, W_PVE, W_PVVip, W_VipSST, D)

# Show raster plot
l1 = Plots.heatmap(v1, title="Layer 1: Somatic Voltage Trace", xaxis="Time (0.1 ms)", yaxis="Nr. of Neuron", c=:balance)
l2 = Plots.heatmap(v2, title="Layer 2: Somatic Voltage Trace", xaxis="Time (0.1 ms)", yaxis="Nr. of Neuron", c=:balance)
display(plot(l1, l2, layout=(2, 1)))
