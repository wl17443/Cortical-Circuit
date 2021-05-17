#=============================================================================================================

Function to simulate interneurons

=============================================================================================================#

function simulateVip(dt, v, u, noise_lvl, Ibg, st_EI, st_II, W_EI, W_II)

    v += dv_vip_dt(v, u, Ibg, st_EI, st_II, W_EI, W_II) * dt
    u += du_vip_dt(v, u) * dt
    Ibg += dIbg_dt(Ibg, dt, noise_lvl) * dt

    newt_ = zeros(length(v))

    # Get spiketiming of each spike
    for i=1:length(v)
        if v[i] >= 50
            newt_[i] = 1
            v[i] = -56
            u[i] = u[i] + 130
        end
    end

    return v, u, Ibg, newt_
end

function simulateSST(dt, v, u, noise_lvl, Ibg, st_EI, st_II, W_EI, W_II)

    v += dv_sst_dt(v, u, Ibg, st_EI, st_II, W_EI, W_II) * dt
    u += du_sst_dt(v, u) * dt
    Ibg += dIbg_dt(Ibg, dt, noise_lvl) * dt

    newt_ = zeros(length(v))

    for i=1:length(v)
        if v[i] >= 40 - 0.1u[i]
            newt_[i] = 1
            v[i] = -53 + 0.04 * u[i]
            u[i] = min(u[i] + 20, 670)
        end
    end

    return v, u, Ibg, newt_
end

function simulatePV(dt, v, u, noise_lvl, Ibg, st_EI, W_EI)

    v += dv_pv_dt(v, u, Ibg, st_EI, W_EI) * dt
    u += du_pv_dt(v, u) * dt
    Ibg += dIbg_dt(Ibg, dt, noise_lvl) * dt

    newt_ = zeros(length(v))

    for i=1:length(v)
        if v[i] >= 25
            newt_[i] = 1
            v[i] = -45
        end
    end

    return v, u, Ibg, newt_
end
