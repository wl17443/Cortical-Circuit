using Plots
using CSV, DataFrames
using Random, Distributions

function find_edges(arr)
    a = 1
    b = 1
    for i=1:length(arr)-1
        arr[i] < arr[i+1] ? a = i+1 : a
        arr[i] > arr[i+1] ? b = i : b
    end
    return a, b
end

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/data/01-04-21-13-41/"
csvfile = readdir(path)[11]

spike_trains = Array(CSV.read("$path$csvfile", DataFrame, header=false))
x = 1:50
means = []
sigmas = []

for i=1:Int(3000/100-1)

y = sum(spike_trains[:, (i-1)*100+1:i*100], dims=2) ./ 100

mean_voltage = mean(y)
min_voltage = minimum(y)
max_voltage = maximum(y)

box = map!(x -> x >= mean_value ? max_value : min_value, y, y)
a, b = find_edges(box)

sigma = a < b ? b - a : length(y) - (a - b)
mu = a + sigma/2 > 50 ? a + sigma/2 - 50 : a + sigma/2

push!(sigmas, sigma)
push!(means, mu)

end

xs = 1:Int(3000/100-1)
display(plot(xs, means, grid=false, yerror=sigmas, ylims=(1, 50)))
