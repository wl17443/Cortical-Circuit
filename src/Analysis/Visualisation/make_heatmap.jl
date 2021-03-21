using Plots
using CSV, DataFrames

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/csv/"
savepath = "C:/Users/Orion/Documents/University/Dissertation/Julia/fig/"
csvfiles = readdir(path)

for csvfile in csvfiles
    v_s = CSV.File("$path$csvfile") |> Tables.matrix
    plt = Plots.heatmap(v_s, c=:balance)
    filename = split(csvfile, ".csv")[1]
    savefig(plt, "$savepath$filename.png")
end
