using LaTeXStrings
using Random, Distributions

function initialise_weights(network_params::Dict, con_params::Dict, noise::Float64)
    Distr_EE = VonMises(0, con_params["k_excexc"])

    W_EE = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pyc"]
        if i != j
            distance = i-j < 0 ? (abs(i-j) <= network_params["nr_pyc"]/2 ? i-j : i-j + network_params["nr_pyc"]) : (i-j >= network_params["nr_pyc"]/2 ? network_params["nr_pyc"]-i+j : i-j)
            # W_EE[i, j] = con_params["s_excexc"] * pdf(Distr_EE, distance * 2*pi/network_params["nr_pyc"]) + rand(Normal(0, con_params["s_excexc"] * con_params["sigma"]))
            W_EE[i, j] = con_params["s_excexc"] * pdf(Distr_EE, distance * 2*pi/network_params["nr_pyc"])
        end
    end

    W_EE += rand(Normal(0, con_params["s_excexc"] * noise), (network_params["nr_pyc"], network_params["nr_pyc"]))

    W_ESST = con_params["s_excsst"] .* (ones(network_params["nr_pyc"], network_params["nr_sst"]) + rand(Normal(0, con_params["s_excsst"] * noise), (network_params["nr_pyc"], network_params["nr_sst"])))
    W_EPV  = con_params["s_excpv"] .* (ones(network_params["nr_pyc"], network_params["nr_pv"]) + rand(Normal(0, con_params["s_excpv"] * noise), (network_params["nr_pyc"], network_params["nr_pv"])))

    W_SSTEd = con_params["s_sstexc"] .* (ones(network_params["nr_sst"], network_params["nr_pyc"]) + rand(Normal(0, con_params["s_sstexc"] * noise), (network_params["nr_sst"], network_params["nr_pyc"])))
    W_PVEs  = con_params["s_pvexc"] .* (ones(network_params["nr_pv"], network_params["nr_pyc"]) + rand(Normal(0, con_params["s_pvexc"] * noise), (network_params["nr_pv"], network_params["nr_pyc"])))

    W_SSTPV = con_params["s_sstpv"] .* (ones(network_params["nr_sst"], network_params["nr_pv"]) + rand(Normal(0, con_params["s_sstpv"] * noise), (network_params["nr_sst"], network_params["nr_pv"])))
    W_PVSST = con_params["s_pvsst"] .* (ones(network_params["nr_pv"], network_params["nr_sst"]) + rand(Normal(0, con_params["s_pvsst"] * noise), (network_params["nr_pv"], network_params["nr_sst"])))

    # W_ESST = con_params["s_excsst"] .* ones(network_params["nr_pyc"], network_params["nr_sst"])
    # W_EPV  = con_params["s_excpv"] .* ones(network_params["nr_pyc"], network_params["nr_pv"])
    #
    # W_SSTEd = con_params["s_sstexc"] .* ones(network_params["nr_sst"], network_params["nr_pyc"])
    # W_PVEs  = con_params["s_pvexc"] .* ones(network_params["nr_pv"], network_params["nr_pyc"])
    #
    # W_SSTPV = con_params["s_sstpv"] .* ones(network_params["nr_sst"], network_params["nr_pv"])
    # W_PVSST = con_params["s_pvsst"] .* ones(network_params["nr_pv"], network_params["nr_sst"])

    return abs.(W_EE), abs.(W_ESST), abs.(W_EPV), abs.(W_SSTEd), abs.(W_PVEs), abs.(W_SSTPV), abs.(W_PVSST)
end
