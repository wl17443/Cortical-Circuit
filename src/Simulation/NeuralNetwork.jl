module NeuralNetwork

include("../Model/Model.jl")
include("../Analysis/BumpCircuitDynamics/BumpDetector.jl")

using .Model: PyramidalNeuron.simulatePyC
using .Model: Interneuron.simulateI
using .BumpDetector: bump_status

include("start_simulation.jl")

export start_simulation

end
