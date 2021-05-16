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
                       "nr_pv"  => 10)

## Simulation timing
t = 10000e-3; dt = 1e-3

## Stimulation
I_inj_s = zeros(network_params["nr_pyc"], Int(t/dt))
stimulation_strength = 5e-6 #nA
stimulation_duration = 50 #ms
add_inj_current!(I_inj_s, stimulation_strength, 1, stimulation_duration, 20)

## Selection
I_inj_d = zeros(network_params["nr_pyc"], Int(t/dt))
selection_strength = 10*stimulation_strength #nA
selection_duration = 5 #ms
add_inj_current!(I_inj_d, selection_strength, 10, selection_duration, 20)

# datetime = Dates.format(now(), "dd-mm-yy-HH-MM-SS")

noise = 0.0
bg_noise = 1

# for s_excpv=14e-9:0.5e-9:18e-9
#
# for bg_noise=[0, 1, 10, 100, 1000, 10000]
#
# sim_log = open("logs/big_noise/s_excpv/$(round(s_excpv*1e9, digits=3))-$bg_noise.txt", "a")
#
# for trial=1:5
#
# write(sim_log, "Trial $trial\n")

# Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "s_excexc" => 3650e-9, # 2900e-9,
                   "s_excsst" => 20e-9,
                   "s_excpv"  => 50e-9,
                   "s_sstexc" => 18e-9,
                   "s_pvexc"  => 80e-9,
                   "s_sstpv"  => 10e-9,
                   "s_pvsst"  => 10e-9,
                   # Concentration parameter kappa
                   "k_excexc" => 60 )

noise_lvl = Dict( "som" => bg_noise,
                  "den" => bg_noise,
                  "inh" => bg_noise)

W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(network_params, con_params, noise)

log = open("results/firint_rates.csv", "w")

v_s, firing_rate = start_simulation(log, t, dt, 500, network_params, noise_lvl, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

firing_rate ./= 0.5

# display(plot(collect(1:20), transpose(firing_rate), legend=false))

# end
#
# close(sim_log)
#
# end
#
# end

# distr = [sum(W_EE)/50 * avg_fr["PyC"], sum(W_ESST)/50 * avg_fr["PyC"], sum(W_EPV)/50 * avg_fr["PyC"], sum(W_SSTE)/5 * avg_fr["SST"], sum(W_PVE)/5 * avg_fr["PV"], sum(W_SSTPV)/5 * avg_fr["SST"], sum(W_PVSST)/5 * avg_fr["PV"]]
# distr = distr ./ sum(distr)
# map!(x -> 100 * x, distr, distr)

# write2csv(v_s, "Voltage-Trace", datetime)

# display(plot(collect(1:200), v_s[20, 1:200]))

Plots.heatmap(v_s, title="Somatic Voltage Trace", xaxis="Time (ms)", yaxis="Nr. of Neuron", c=:balance)
# Plots.heatmap(W_EE)
