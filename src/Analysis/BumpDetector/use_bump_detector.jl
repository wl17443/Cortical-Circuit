include("bump_detector_funcs.jl")

using Plots
using CSV, DataFrames

df = DataFrame(CSV.File("C:/Users/Orion/Documents/University/Dissertation/Julia/data/Voltage-Trace/03-05-21-15-20-06-voltage.csv", header=false))
vol = Matrix(df)

# activity = sum(vol[:, 1:1000], dims=2) ./ 1000
# meanv = round(mean(activity), digits=4)
# minv = round(minimum(activity), digits=4)
# maxv = round(maximum(activity), digits=4)
sigma, mu = bump_status(vol[:, 1:1000], 1000, 50)

# xs = collect(1:50)
# plt  = plot(xs, activity, legend=false, margin=1Plots.mm)
# plot!([meanv, maxv, minv], seriestype="hline", yticks = ([meanv, maxv, minv],["$meanv","$maxv", "$minv"]), line=:dash, label=false)
# title!("Average Voltage of all PyCs across 1000ms of Simulation")
# xlabel!("Neuron Number")
# ylabel!("Voltage (mV)")

# savefig(plt, "C:/Users/Orion/Documents/University/Dissertation/Figs-Results/voltage-breakdown-for-eval-section.png")

# anim = @animate for t in 1:50
#     plot(xs, vol[:, t], title="Time: $t ms", legend=false)
# end

# gif(anim, "test.gif", fps=3)
