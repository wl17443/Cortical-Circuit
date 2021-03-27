#=============================================================================================================

Function to initialise injected current

=============================================================================================================#

function initialise_inj_current(nr_neurons::Int, steps::Int, amount::Float64, start::Int, duration::Int)
    inj_current = zeros(nr_neurons, steps)
    inj_current[10, start:(start+duration)] .= amount

    return inj_current
end

function add_inj_current!(inj_current, amount, start, duration, at_neuron)
    inj_current[at_neuron, start:(start+duration)] .= amount
end
