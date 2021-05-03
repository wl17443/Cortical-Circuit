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
                       "nr_pyc" => 50,
                       "nr_sst" => 5,
                       "nr_pv"  => 5)

## Simulation timing
t = 6000e-3; dt = 1e-3

## Stimulation
I_inj_s = zeros(network_params["nr_pyc"], Int(t/dt))
stimulation_strength = 5e-6 #nA
stimulation_duration = 50 #ms
add_inj_current!(I_inj_s, stimulation_strength, 100, stimulation_duration, 10)

## Selection
I_inj_d = zeros(network_params["nr_pyc"], Int(t/dt))
selection_strength = 10*stimulation_strength #nA
selection_duration = 5 #ms
add_inj_current!(I_inj_d, selection_strength, 110, selection_duration, 10)

# datetime = Dates.format(now(), "dd-mm-yy-HH-MM-SS")

noise = 0.0

for k_excexc=0.5:0.5:20

sim_log = open("logs/parameter_sweep/k_excexc/$(round(k_excexc, digits=3)).txt", "a")

for trial=1:5

write(sim_log, "Trial $trial\n")

# Connectivity
con_params = Dict( # Post synaptic current amplitude
                   "s_excexc" => 580e-9,
                   "s_excsst" => 18e-9,
                   "s_excpv"  => 18e-9,
                   "s_sstexc" => 15e-9,
                   "s_pvexc"  => 65e-9,
                   "s_sstpv"  => 10e-9,
                   "s_pvsst"  => 10e-9,
                   # Concentration parameter kappa
                   "k_excexc" => k_excexc )

noise_lvl = Dict( "som" => 1,
                  "den" => 1,
                  "inh" => 1)

W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST = initialise_weights(network_params, con_params, noise)

_ = start_simulation(sim_log, t, dt, 600, network_params, noise_lvl, I_inj_s, I_inj_d, W_EE, W_ESST, W_EPV, W_SSTE, W_PVE, W_SSTPV, W_PVSST)

end

close(sim_log)

end

# distr = [sum(W_EE)/50 * avg_fr["PyC"], sum(W_ESST)/50 * avg_fr["PyC"], sum(W_EPV)/50 * avg_fr["PyC"], sum(W_SSTE)/5 * avg_fr["SST"], sum(W_PVE)/5 * avg_fr["PV"], sum(W_SSTPV)/5 * avg_fr["SST"], sum(W_PVSST)/5 * avg_fr["PV"]]
# distr = distr ./ sum(distr)
# map!(x -> 100 * x, distr, distr)

# write2csv(v_s, "Voltage-Trace", datetime)
# show_heatmap(v_s, "v_s")
# save_heatmap(v_s, "Noise vs. drift: base+a_pvexc=$(con_params["a_pvexc"]), noise_lvl=$noise_lvl", "noise_pvexc2", "$(con_params["a_pvexc"])-$noise_lvl-$trial")

# if error < 10
#     writedlm(results, con_params)
#
#     println("Synaptic Weights")
#     println(con_params)
#     println("Error")
#     println(error)
# end

# end

# close(results)
