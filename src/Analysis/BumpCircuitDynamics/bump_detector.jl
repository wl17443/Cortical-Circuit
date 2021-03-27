using Plots
using CSV, DataFrames
using Random, Distributions

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/csv/two_activities/"
csvfile = readdir(path)[7]

v_s = CSV.File("$path$csvfile") |> Tables.matrix
mean = zeros(2000)
for step in 1:2000
    norm = fit(Normal, v_s[:, step])
    mean[step], _ = params(norm)
end
display(Plots.heatmap(v_s, c=:balance))
display(Plots.plot(1:2000, mean))
