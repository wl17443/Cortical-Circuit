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

## Stimulation
I_inj_s = zeros(network_params["nr_pyc"], Int(t/dt))
stimulation_strength = 5e-6 #nA
stimulation_duration = 50 #ms
add_inj_current!(I_inj_s, stimulation_strength, 1, stimulation_duration, 10)

## Selection
I_inj_d = zeros(network_params["nr_pyc"], Int(t/dt))
selection_strength = 10*stimulation_strength #nA
selection_duration = 5 #ms
add_inj_current!(I_inj_d, selection_strength, 10, selection_duration, 10)

# datetime = Dates.format(now(), "dd-mm-yy-HH-MM-SS")

# noise = [1; collect(10:10:100)]
# distr = zeros(length(kappa), 7)

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

# for nsom in [0; 1; collect(1000:1000:10000)], nden in [0; 1; collect(1000:1000:10000)], ninh in [0; 1; collect(1000:1000:10000)]

noise_lvl = Dict( "som" => 0,
                  "den" => 0,
                  "inh" => 0)

sim_log = open("logs/distr_test/test4.txt", "a")

# writedlm(sim_log, noise_lvl)
W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(network_params, con_params)

# write(sim_log, "Trial $trial\n")

v_s, distr, firing_rate = start_simulation(sim_log, t, dt, 500, network_params, noise_lvl, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

close(sim_log)

# distr = [sum(W_EE)/50 * avg_fr["PyC"], sum(W_ESST)/50 * avg_fr["PyC"], sum(W_EPV)/50 * avg_fr["PyC"], sum(W_SSTE)/5 * avg_fr["SST"], sum(W_PVE)/5 * avg_fr["PV"], sum(W_SSTPV)/5 * avg_fr["SST"], sum(W_PVSST)/5 * avg_fr["PV"]]
# distr = distr ./ sum(distr)
# map!(x -> 100 * x, distr, distr)

results = open("results.csv", "w")
writedlm(results, firing_rate)
writedlm(results, distr)
close(results)

# write2csv(v_s, "big_noise_base_params", "")
show_heatmap(v_s, "v_s")
# save_heatmap(v_s, "Noise vs. drift: base+a_pvexc=$(con_params["a_pvexc"]), noise_lvl=$noise_lvl", "noise_pvexc2", "$(con_params["a_pvexc"])-$noise_lvl-$trial")

# end
