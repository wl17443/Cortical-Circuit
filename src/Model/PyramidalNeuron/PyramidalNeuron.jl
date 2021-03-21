module PyramidalNeuron
    include("SomaticCompartment.jl")
    include("DendriticCompartment.jl")
    include("simulatePyc.jl")

    export simulatePyc
end
