#=============================================================================================================

Function to simulate interneuron. 

=============================================================================================================#

function simulateVip(dt, v, u, st_E1I, st_E2I, st_II, W_EI, W_II, D)

    v += dv_vip_dt(v, u, st_E1I, st_E2I, st_II, W_EI, W_II, D, dt) * dt
    u += du_vip_dt(v, u) * dt

    newt_ = zeros(length(v))

    # Get spiketiming of each spike
    for i=1:length(v)
        if v[i] >= 50
            newt_[i] = 1
            v[i] = -56
            u[i] = u[i] + 130
        end
    end

    return v, u, newt_
end

function simulateSST(dt, v, u, st_E1I, st_E2I, st_II, W_EI, W_II, D)

    v += dv_sst_dt(v, u, st_E1I, st_E2I, st_II, W_EI, W_II, D, dt) * dt
    u += du_sst_dt(v, u) * dt

    newt_ = zeros(length(v))

    for i=1:length(v)
        if v[i] >= 40 - 0.1u[i]
            newt_[i] = 1
            v[i] = -53 + 0.04 * u[i]
            u[i] = min(u[i] + 20, 670)
        end
    end

    return v, u, newt_
end

function simulatePV(dt, v, u, st_E1I, st_E2I, W_EI, D)

    v += dv_pv_dt(v, u, st_E1I, st_E2I, W_EI, D, dt) * dt
    u += du_pv_dt(v, u) * dt

    newt_ = zeros(length(v))

    for i=1:length(v)
        if v[i] >= 25
            newt_[i] = 1
            v[i] = -45
        end
    end

    return v, u, newt_
end
