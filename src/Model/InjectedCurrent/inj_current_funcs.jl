#=============================================================================================================

Functions to initialise injected current

=============================================================================================================#

function initialise_inj_current(stimulus_params::Dict, network_params::Dict, steps::Int)
    inj_current = zeros(network_params["nr_pyc"], steps)

    if stimulus_params["type"] == "rotation"

        inj_current[at_neuron, start:(start+duration)] .= amount
    end

    return inj_current
end

function add_inj_current!(inj_current::Array{Float64,2}, amount, start, duration, at_neuron)
    inj_current[at_neuron, start:(start+duration)] .= amount
end
