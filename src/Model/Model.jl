module Model
    include("Connectivity/Connectivity.jl")
    include("Interneuron/Interneuron.jl")
    include("PyramidalNeuron/PyramidalNeuron.jl")
    include("InjectedCurrent/InjectedCurrent.jl")

    export Connectivity
    export Interneuron
    export PyramidalNeuron
    export InjectedCurrent
end
