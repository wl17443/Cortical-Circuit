using CSV
using DataFrames
using Plots

df = CSV.read("C:/Users/Orion/Documents/University/Dissertation/Julia/results.csv", DataFrame)
results = Matrix(results)
filtered_som = results[(results[:, 1] .== 1), 4]
filtered_den = results[(results[:, 2] .== 1), 4]
filtered_inh = results[(results[:, 3] .== 1), 4]

ranges = [0; 1; collect(1000:1000:10000)]
len = length(ranges)

reshaped_som = reshape(filtered_som, (len, len))
reshaped_den = reshape(filtered_den, (len, len))
reshaped_inh = reshape(filtered_inh, (len, len))

# print(size(filtered_som))
Plots.heatmap(ranges, ranges, reshaped_inh)
