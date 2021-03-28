using Plots
using CSV, DataFrames

path = "/home/anhelka/Documents/Cortical-Circuit/data/kappa_globalvlocal/"
savepath = "/home/anhelka/Documents/Cortical-Circuit/figs/kappa_globalvlocal/"
csvfiles = readdir(path)

for csvfile in csvfiles
    v_s = CSV.File("$path$csvfile") |> Tables.matrix
    plt = Plots.heatmap(v_s, c=:balance)
    filename = split(csvfile, ".csv")[1]
    savefig(plt, "$savepath$filename.png")
end
