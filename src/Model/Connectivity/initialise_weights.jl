using LaTeXStrings
using Random, Distributions

function initialise_weights(network_params::Dict, con_params::Dict)
    Distr_EE = VonMises(0, con_params["k_excexc"])

    W_EE = zeros(network_params["nr_pyc"], network_params["nr_pyc"])

    for i=1:network_params["nr_pyc"], j=1:network_params["nr_pyc"]
        if i != j
            distance = i-j < 0 ? (abs(i-j) <= network_params["nr_pyc"]/2 ? i-j : i-j + network_params["nr_pyc"]) : (i-j >= network_params["nr_pyc"]/2 ? network_params["nr_pyc"]-i+j : i-j)
            W_EE[i, j] = con_params["a_excexc"] * pdf(Distr_EE, distance * 2*pi/network_params["nr_pyc"])
        end
    end

    # W_ESST = con_params["a_excsst"] * ones(network_params["nr_pyc"], network_params["nr_sst"]) ./ (network_params["nr_pyc"] * network_params["nr_sst"])
    # W_EPV  = con_params["a_excpv"] * ones(network_params["nr_pyc"], network_params["nr_pv"]) ./ (network_params["nr_pyc"] * network_params["nr_pv"])
    #
    # W_SSTEd = con_params["a_sstexc"] * ones(network_params["nr_sst"], network_params["nr_pyc"]) ./ (network_params["nr_sst"] * network_params["nr_pyc"])
    # W_PVEs  = con_params["a_pvexc"] * ones(network_params["nr_pv"], network_params["nr_pyc"]) ./ (network_params["nr_pv"] * network_params["nr_pyc"])
    #
    # W_SSTPV = con_params["a_sstpv"] * ones(network_params["nr_sst"], network_params["nr_pv"]) ./ (network_params["nr_sst"] * network_params["nr_pv"])
    # W_PVSST  = con_params["a_pvsst"] * ones(network_params["nr_pv"], network_params["nr_sst"]) ./ (network_params["nr_pv"] * network_params["nr_sst"])

    W_ESST = rand(Normal(con_params["a_excsst"], 0), (network_params["nr_pyc"], network_params["nr_sst"])) ./ (network_params["nr_pyc"] * network_params["nr_sst"])
    W_EPV  = rand(Normal(con_params["a_excpv"], 0), (network_params["nr_pyc"], network_params["nr_pv"])) ./ (network_params["nr_pyc"] * network_params["nr_pv"])

    W_SSTEd = rand(Normal(con_params["a_sstexc"], 0), (network_params["nr_sst"], network_params["nr_pyc"])) ./ (network_params["nr_sst"] * network_params["nr_pyc"])
    W_PVEs  = rand(Normal(con_params["a_pvexc"], 0), (network_params["nr_pv"], network_params["nr_pyc"])) ./ (network_params["nr_pv"] * network_params["nr_pyc"])

    W_SSTPV = rand(Normal(con_params["a_sstpv"], 0), (network_params["nr_sst"], network_params["nr_pv"])) ./ (network_params["nr_sst"] * network_params["nr_pv"])
    W_PVSST = rand(Normal(con_params["a_pvsst"], 0), (network_params["nr_pv"], network_params["nr_sst"])) ./ (network_params["nr_pv"] * network_params["nr_sst"])

    return W_EE, W_ESST, W_EPV, W_SSTEd, W_PVEs, W_SSTPV, W_PVSST
end
