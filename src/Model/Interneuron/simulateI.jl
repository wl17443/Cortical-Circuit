#=============================================================================================================

Function to simulate interneurons

=============================================================================================================#

v_thr = -44.5e-3; v_rest = -55.5e-3

function simulateI(t, dt, v, noise_lvl, Ibg, tspike, st_EI, st_II, W_EI, W_II)
    v += dv_dt(v, Ibg, st_EI, st_II, W_EI, W_II) .* dt
    Ibg += dIbg_dt(Ibg, dt, noise_lvl) .* dt

    # Get spike train for this timestep
    newt_ = map(x -> x >= v_thr ? 1 : 0, v)
    # Get spiketiming of each spike
    for i=1:length(tspike)
        if t * dt - tspike[i] <= 3e-3
            v[i] = v_rest
        end
        if newt_[i] == 1
            tspike[i] = t * dt
            v[i] = v_rest
        end
    end

    return v, Ibg, newt_, tspike
end
