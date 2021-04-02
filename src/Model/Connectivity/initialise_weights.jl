using Random, Distributions

function initialise_weights(nr_pyc, nr_sst, nr_pv, con_params::Dict)
    ## Connectivity matrix that makes up the ring attractor network of PyC
    Distr_EE = VonMises(0, con_params["kappa_excexc"])

    W_EE = zeros(nr_pyc, nr_pyc)

    for i=1:nr_pyc, j=1:nr_pyc
        if i != j
            distance = i-j < 0 ? (abs(i-j) <= nr_pyc/2 ? i-j : i-j + nr_pyc) : (i-j >= nr_pyc/2 ? nr_pyc-i+j : i-j)
            W_EE[i, j] = con_params["alpha_excexc"] * pdf(Distr_EE, distance * 2*pi/nr_pyc)
        end
    end

    ## Synaptic Weights (Fixed)
    ## E->I{SST,PV}
    # Distr_ESST = VonMises(0, con_params["kappa_excinh"])
    # Distr_EPV  = VonMises(0, con_params["kappa_excinh"])

    W_ESST = con_params["alpha_excinh"] * ones(nr_pyc, nr_sst) ./ (nr_pyc * nr_sst)
    W_EPV  = con_params["alpha_excinh"] * ones(nr_pyc, nr_pv) ./ (nr_pyc * nr_pv)

    # for i=1:nr_pyc, j=1:nr_sst
    #     distance = i-j < 0 ? (abs(i-j) <= nr_sst/2 ? i-j : i-j + nr_sst) : (i-j >= nr_sst/2 ? nr_sst-i+j : i-j)
    #     W_ESST[i, j] = con_params["alpha_excinh"] * pdf(Distr_ESST, distance * 2*pi/nr_sst)
    # end
    #
    # for i=1:nr_pyc, j=1:nr_pv
    #     distance = i-j < 0 ? (abs(i-j) <= nr_pv/2 ? i-j : i-j + nr_pv) : (i-j >= nr_pv/2 ? nr_pv-i+j : i-j)
    #     W_EPV[i, j] = con_params["alpha_excinh"] * pdf(Distr_EPV, distance * 2*pi/nr_pv)
    # end

    ## I{SST,PV}->E{s,d}
    # Distr_SSTEd = VonMises(0, con_params["kappa_inhexc"])
    # Distr_PVEs  = VonMises(0, con_params["kappa_inhexc"])

    W_SSTEd = con_params["alpha_inhexc"] * ones(nr_sst, nr_pyc) ./ (nr_sst * nr_pyc)
    W_PVEs  = con_params["alpha_inhexc"] * ones(nr_pv, nr_pyc) ./ (nr_pv * nr_pyc)

    # for i=1:nr_sst, j=1:nr_pyc
    #     distance = i-j < 0 ? (abs(i-j) <= nr_pyc/2 ? i-j : i-j + nr_pyc) : (i-j >= nr_pyc/2 ? nr_pyc-i+j : i-j)
    #     W_SSTEd[i, j] = con_params["alpha_inhexc"] * pdf(Distr_SSTEd, distance * 2*pi/nr_pyc)
    # end
    #
    # for i=1:nr_pv, j=1:nr_pyc
    #     distance = i-j < 0 ? (abs(i-j) <= nr_pyc/2 ? i-j : i-j + nr_pyc) : (i-j >= nr_pyc/2 ? nr_pyc-i+j : i-j)
    #     W_PVEs[i, j] = con_params["alpha_inhexc"] * pdf(Distr_PVEs, distance * 2*pi/nr_pyc)
    # end

    ## I{SST,PV}->I{SST,PV}
    # Distr_SSTPV = VonMises(0, con_params["kappa_inhinh"])
    # Distr_PVSST  = VonMises(0, con_params["kappa_inhinh"])

    W_SSTPV = con_params["alpha_inhinh"] * ones(nr_sst, nr_pv) ./ (nr_sst * nr_pv)
    W_PVSST  = con_params["alpha_inhinh"] * ones(nr_pv, nr_sst) ./ (nr_pv * nr_sst)

    # for i=1:nr_sst, j=1:nr_pv
    #     distance = i-j < 0 ? (abs(i-j) <= nr_pv/2 ? i-j : i-j + nr_pv) : (i-j >= nr_pv/2 ? nr_pv-i+j : i-j)
    #     W_SSTPV[i, j] = con_params["alpha_inhinh"] * pdf(Distr_SSTPV, distance * 2*pi/nr_pv)
    # end
    #
    # for i=1:nr_pv, j=1:nr_sst
    #     distance = i-j < 0 ? (abs(i-j) <= nr_sst/2 ? i-j : i-j + nr_sst) : (i-j >= nr_sst/2 ? nr_sst-i+j : i-j)
    #     W_PVSST[i, j] = con_params["alpha_inhinh"] * pdf(Distr_PVSST, distance * 2*pi/nr_sst)
    # end

    return W_EE, W_ESST, W_EPV, W_SSTEd, W_PVEs, W_SSTPV, W_PVSST
end
