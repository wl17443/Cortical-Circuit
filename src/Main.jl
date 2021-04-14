include("Simulation/NeuralNetwork.jl")
include("Analysis/WriteData/WriteData.jl")
include("Analysis/Visualisation/Visualise.jl")

using CSV
using Plots
using Dates
using DelimitedFiles
using .WriteData: write2csv
using .Visualise: save_heatmap, show_heatmap
using .NeuralNetwork: Model.Connectivity.initialise_weights
using .NeuralNetwork: start_simulation
using .NeuralNetwork: Model.InjectedCurrent.initialise_inj_current, Model.InjectedCurrent.add_inj_current!

network_params = Dict( # Network parameters - size
                       "nr_pyc" => 50,
                       "nr_sst" => 5,
                       "nr_pv"  => 5)

## Simulation timing
t = 10000e-3; dt = 1e-3

# stimulus_params = Dict( "type" => "rotation",
#     	                "start_neuron" => 1,
#                         "rate" => 1, # neuron per ms
#                         "direction" => "left",
#                         "strength" => 5e-6)

# stimulus_params = Dict( "type" => "single",
#                         "start_neuron" => 1,
#                         "rate" => 1, # neuron per ms
#                         "direction" => "left",
#                         "strength" => 5e-6)

## Stimulation
I_inj_s = zeros(network_params["nr_pyc"], Int(t/dt))
stimulation_strength = 5e-6 #nA
stimulation_duration = 50 #ms
add_inj_current!(I_inj_s, stimulation_strength, 1, stimulation_duration, 10)
add_inj_current!(I_inj_s, stimulation_strength, 2000, stimulation_duration, 20)

## Selection
I_inj_d = zeros(network_params["nr_pyc"], Int(t/dt))
selection_strength = 10*stimulation_strength #nA
selection_duration = 5 #ms
add_inj_current!(I_inj_d, selection_strength, 10, selection_duration, 10)
add_inj_current!(I_inj_d, selection_strength, 2010, selection_duration, 20)

datetime = Dates.format(now(), "dd-mm-yy-HH-MM-SS")

log = open("logs/13-04-2021-16-58/$datetime.txt", "a")

# Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "a_excexc" => 6e-7,
                   "a_excsst" => 45e-7,
                   "a_excpv" => 45e-7,
                   "a_sstexc" => 10e-7,
                   "a_pvexc" => 200e-7,
                   "a_sstpv" => 0.1e-7,
                   "a_pvsst" => 0.1e-7,
                   # Concentration parameter kappa
                   "k_excexc" => 13.5)

writedlm(log, con_params)
W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(network_params, con_params)

write(log, datetime * "\n")

v_s = start_simulation(log, t, dt, 1000, network_params, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

write2csv(v_s, "13-04-2021-16-58", datetime)
show_heatmap(v_s, "Opposing Intermittent Stimuli: base parameters, a_ssteexc=$(con_params["a_sstexc"])")
save_heatmap(v_s, "Opposing Intermittent Stimuli: base parameters, a_ssteexc=$(con_params["a_sstexc"])", "13-04-2021-16-58", datetime)

close(log)
