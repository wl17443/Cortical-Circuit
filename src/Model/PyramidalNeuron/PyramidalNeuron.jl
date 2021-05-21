module PyramidalNeuron
    include("somatic_funcs.jl")
    include("dendritic_funcs.jl")
    include("simulatePyC.jl")

    export simulatePyC1, simulatePyC2
end
