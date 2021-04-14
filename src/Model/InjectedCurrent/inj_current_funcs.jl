#=============================================================================================================

Functions to initialise injected current

=============================================================================================================#

function initialise_inj_current(stimulus_params::Dict, network_params::Dict, steps::Int)
    inj_current = zeros(network_params["nr_pyc"], steps)

    # stimulus_params = Dict( "type" => "rotation",
    #     	                "start_neuron" => 1,
    #                         "rate" => 1, # neuron per ms
    #                         "direction" => "left",
    #                         "strength" => 5e-6)

    if stimulus_params["type"] == "rotation"
        
        inj_current[at_neuron, start:(start+duration)] .= amount
    end

    return inj_current
end

function add_inj_current!(inj_current::Array{Float64,2}, amount, start, duration, at_neuron)
    inj_current[at_neuron, start:(start+duration)] .= amount
end
