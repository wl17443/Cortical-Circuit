using Plots
using CSV, DataFrames
using Random, Distributions
using KernelDensity

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/data/kappa_globalvlocal_spiketrains/"
csvfile = readdir(path)[1]

spike_trains = Array(CSV.read("$path$csvfile", DataFrame, header=false))
x = spike_trains[:, 2100]
indices = 1:50
sum_indices = sum(indices)

data_plot = x .* indices
gauss_data = fit(Normal, data_plot)

# estimation = kde(x, kernel=Normal)
# x_pdf = pdf(estimation, x)

print(params(gauss_data))
