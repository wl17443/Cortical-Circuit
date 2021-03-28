using Plots
using CSV, DataFrames
using Random, Distributions
using KernelDensity

path = "/home/anhelka/Documents/Cortical-Circuit/data/kappa_globalvlocal_spiketrains/"
csvfile = readdir(path)[9]

spike_trains = Array(CSV.read("$path$csvfile", DataFrame, header=false))
x = spike_trains[:, 50]
estimation = kde(x, kernel=Normal)
x_pdf = pdf(estimation, x)

display(Plots.plot(1:50, x_pdf))
