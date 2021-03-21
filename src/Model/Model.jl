module Model    
    include("Connectivity/Connectivity.jl")
    include("Interneuron/Interneuron.jl")
    include("PyramidalNeuron/Pyramidalneuron.jl")

    export Connectivity
    export Interneuron
    export PyramidalNeuron
end
