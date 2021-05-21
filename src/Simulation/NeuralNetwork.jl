module NeuralNetwork

include("../Model/Model.jl")
include("../Analysis/BumpDetector/BumpDetector.jl")

using .Model: PyramidalNeuron.simulatePyC1, PyramidalNeuron.simulatePyC2
using .Model: Interneuron.simulateSST, Interneuron.simulatePV, Interneuron.simulateVip
using .BumpDetector: bump_status

include("start_simulation.jl")

export start_simulation

end
