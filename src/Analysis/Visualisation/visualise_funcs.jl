#=============================================================================================================

Functions used for visualising results. 

=============================================================================================================#

using Plots
using Dates

default(titlefont=(10, "sans-serif"), legendfontsize=18)

function show_heatmap(arr::Array{Float64,2}, title)
    plt = Plots.heatmap(arr, c=:balance, title=title)
    display(plt)
end

function csv2heatmap(input_folder::String, save_folder::String)
    path = "C:/Users/Orion/Documents/University/Dissertation/Julia/data/$input_folder/"
    savepath = "C:/Users/Orion/Documents/University/Dissertation/Julia/figs/$save_folder/"
    csvfiles = readdir(path)

    for csvfile in csvfiles
        v_s = Array(CSV.read("$path$csvfile", DataFrame, header=false))
        plt = Plots.heatmap(v_s, c=:balance)
        filename = split(csvfile, ".csv")[1]
        savefig(plt, "$savepath$filename.png")
    end
end

function save_heatmap(arr::Array{Float64,2}, title::String, folder_name::String, datetime::String)
    plt = Plots.heatmap(arr, c=:balance, title=title, xlabel="Time (ms)", ylabel="Neuron Nr.")

    savepath = "C:/Users/Orion/Documents/University/Dissertation/Julia/data/$folder_name/"
    savefig(plt, "$savepath$datetime.png")
end
