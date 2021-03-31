using Plots
using CSV, DataFrames

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/data/volt-kappa_excexc/"
savepath = "C:/Users/Orion/Documents/University/Dissertation/Julia/figs/volt-kappa_excexc/"
csvfiles = readdir(path)

for csvfile in csvfiles
    v_s = Array(CSV.read("$path$csvfile", DataFrame, header=false))
    plt = Plots.heatmap(v_s, c=:balance)
    filename = split(csvfile, ".csv")[1]
    savefig(plt, "$savepath$filename.png")
end
