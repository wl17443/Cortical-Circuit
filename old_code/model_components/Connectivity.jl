module Connectivity

using ..Units 
using Random, Distributions
using Plots

export intialise_weights

function initialise_weights(nr_pyc, nr_sst, nr_pv, alpha_exc, alpha_inh)

## Connectivity matrix that makes up the ring attractor network of PyC
Distr_EE = VonMises(0, 1/0.168^2)

W_EE = zeros(nr_pyc, nr_pyc)

for i=1:nr_pyc, j=1:nr_pyc
    if i != j
        distance = i-j < 0 ? (abs(i-j) <= nr_pyc/2 ? i-j : i-j + nr_pyc) : (i-j >= nr_pyc/2 ? nr_pyc-i+j : i-j)
        W_EE[i, j] = alpha_exc * pdf(Distr_EE, distance * 2*pi/nr_pyc)
    end
end

## Synaptic Weights (Fixed)
## E->I{SST,PV}
Distr_ESST = VonMises(0, 1/0.244^2)
Distr_EPV  = VonMises(0, 1/0.575^2)

W_ESST = zeros(nr_pyc, nr_sst)
W_EPV  = zeros(nr_pyc, nr_pv)

for i=1:nr_pyc, j=1:nr_sst
    if i != j
        distance = i-j < 0 ? (abs(i-j) <= nr_sst/2 ? i-j : i-j + nr_sst) : (i-j >= nr_sst/2 ? nr_sst-i+j : i-j)
        W_ESST[i, j] = alpha_exc * pdf(Distr_ESST, distance * 2*pi/nr_sst)
    end
end

for i=1:nr_pyc, j=1:nr_pv
    if i != j
        distance = i-j < 0 ? (abs(i-j) <= nr_pv/2 ? i-j : i-j + nr_pv) : (i-j >= nr_pv/2 ? nr_pv-i+j : i-j)
        W_EPV[i, j] = alpha_exc * pdf(Distr_EPV, distance * 2*pi/nr_pv)
    end
end

## I{SST,PV}->E{s,d}
Distr_SSTEd = VonMises(0, 1/0.465^2)
Distr_PVEs  = VonMises(0, 1/0.600^2)

W_SSTEd = zeros(nr_sst, nr_pyc)
W_PVEs  = zeros(nr_pv, nr_pyc)

for i=1:nr_sst, j=1:nr_pyc
    if i != j
        distance = i-j < 0 ? (abs(i-j) <= nr_pyc/2 ? i-j : i-j + nr_pyc) : (i-j >= nr_pyc/2 ? nr_pyc-i+j : i-j)
        W_SSTEd[i, j] = alpha_inh * pdf(Distr_SSTEd, distance * 2*pi/nr_pyc)
    end
end

for i=1:nr_pv, j=1:nr_pyc
    if i != j
        distance = i-j < 0 ? (abs(i-j) <= nr_pyc/2 ? i-j : i-j + nr_pyc) : (i-j >= nr_pyc/2 ? nr_pyc-i+j : i-j)
        W_PVEs[i, j] = alpha_inh * pdf(Distr_PVEs, distance * 2*pi/nr_pyc)
    end
end

## I{SST,PV}->I{SST,PV}
Distr_SSTPV = VonMises(0, 1/0.379^2)
Distr_PVSST  = VonMises(0, 1/0.241^2)

W_SSTPV = zeros(nr_sst, nr_pyc)
W_PVSST  = zeros(nr_pv, nr_pyc)

for i=1:nr_sst, j=1:nr_pv
    if i != j
        distance = i-j < 0 ? (abs(i-j) <= nr_pv/2 ? i-j : i-j + nr_pv) : (i-j >= nr_pv/2 ? nr_pv-i+j : i-j)
        W_SSTPV[i, j] = alpha_inh * pdf(Distr_SSTPV, distance * 2*pi/nr_pv)
    end
end

for i=1:nr_pv, j=1:nr_sst
    if i != j
        distance = i-j < 0 ? (abs(i-j) <= nr_sst/2 ? i-j : i-j + nr_sst) : (i-j >= nr_sst/2 ? nr_sst-i+j : i-j)
        W_PVSST[i, j] = alpha_inh * pdf(Distr_PVSST, distance * 2*pi/nr_sst)
    end
end

return W_EE, W_ESST, W_EPV, W_SSTEd, W_PVEs, W_SSTPV, W_PVSST

end

end
