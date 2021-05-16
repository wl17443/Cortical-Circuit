using Statistics
using Plots

function find_edges(arr)
    a = 1
    b = 1
    for i=1:length(arr)-1
        arr[i] < arr[i+1] ? a = i+1 : a
        arr[i] > arr[i+1] ? b = i : b
    end
    return a, b
end

function bump_status(spike_trains, slice_size, ring_size)
    activity = sum(spike_trains, dims=2) ./ slice_size

    mean_voltage = mean(activity)
    min_voltage = minimum(activity)
    max_voltage = maximum(activity)

    box = map!(x -> x >= mean_voltage ? max_voltage : min_voltage, activity, activity)

    # plt = plot(collect(1:50), box, legend=false)
    # title!("Boxed Result")
    # xlabel!("Neuron Number")
    # ylabel!("Voltage (mV)")
    # savefig(plt, "C:/Users/Orion/Documents/University/Dissertation/Figs-Results/voltage-box-for-eval-section.png")

    a, b = find_edges(box)

    sigma = a < b ? b - a : ring_size - (a - b)
    mu = a + sigma/2 > ring_size ? a + sigma/2 - ring_size : a + sigma/2

    return sigma, mu
end
