module NeuralNetwork

include("../Model/Model.jl")

using .Model: PyramidalNeuron.simulatePyC
using .Model: Interneuron.simulateI

include("start_simulation.jl")

export start_simulation

end
