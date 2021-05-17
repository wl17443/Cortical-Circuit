using LaTeXStrings
using LinearAlgebra
using Random, Distributions

function initialise_weights(network_params::Dict, con_params::Dict, noise::Float64)
    Distr_EE = VonMises(0, con_params["k_excexc"])
    weighted = pdf.(Distr_EE, collect(0:(2pi/network_params["nr_pyc"]):pi))
    weighted /= sum(weighted)

    W_EE = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pyc"]
        if i != j
            distance = abs(i-j) > network_params["nr_pyc"]/2 ? network_params["nr_pyc"] - abs(i-j) : abs(i-j)
            W_EE[i, j] = con_params["s_excexc"] * weighted[distance]
        end
    end

    W_EE += rand(Normal(0, con_params["s_excexc"] * noise), (network_params["nr_pyc"], network_params["nr_pyc"]))

    W_ESST = con_params["s_excsst"] .* (ones(network_params["nr_pyc"], network_params["nr_sst"]) + rand(Normal(0, con_params["s_excsst"] * noise), (network_params["nr_pyc"], network_params["nr_sst"])))
    W_EPV  = con_params["s_excpv"] .* (ones(network_params["nr_pyc"], network_params["nr_pv"]) + rand(Normal(0, con_params["s_excpv"] * noise), (network_params["nr_pyc"], network_params["nr_pv"])))
    W_EVip  = con_params["s_excvip"] .* (ones(network_params["nr_pyc"], network_params["nr_vip"]) + rand(Normal(0, con_params["s_excvip"] * noise), (network_params["nr_pyc"], network_params["nr_vip"])))

    W_SSTEd = con_params["s_sstexc"] .* (ones(network_params["nr_sst"], network_params["nr_pyc"]) + rand(Normal(0, con_params["s_sstexc"] * noise), (network_params["nr_sst"], network_params["nr_pyc"])))
    W_PVEs  = con_params["s_pvexc"] .* (ones(network_params["nr_pv"], network_params["nr_pyc"]) + rand(Normal(0, con_params["s_pvexc"] * noise), (network_params["nr_pv"], network_params["nr_pyc"])))

    W_PVVip = con_params["s_pvvip"] .* (ones(network_params["nr_pv"], network_params["nr_vip"]) + rand(Normal(0, con_params["s_pvvip"] * noise), (network_params["nr_pv"], network_params["nr_vip"])))
    W_VipSST = con_params["s_vipsst"] .* (ones(network_params["nr_vip"], network_params["nr_sst"]) + rand(Normal(0, con_params["s_vipsst"] * noise), (network_params["nr_vip"], network_params["nr_sst"])))

    return abs.(W_EE), abs.(W_ESST), abs.(W_EPV), abs.(W_EVip), abs.(W_SSTEd), abs.(W_PVEs), abs.(W_PVVip), abs.(W_VipSST)
end
