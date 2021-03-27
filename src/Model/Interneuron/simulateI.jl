#=============================================================================================================

Function to simulate interneurons 

=============================================================================================================#

function simulateI(t, dt, v, Ibg, tspike, st_EI, st_II, W_EI, W_II)
    v += dv_dt(v, Ibg, st_EI, st_II, W_EI, W_II) .* dt
    Ibg += dIbg_dt(Ibg, dt) .* dt

    # Get spike train for this timestep
    newt_ = map(x -> x >= v_thr ? 1 : 0, v)
    # Get spiketiming of each spike
    for i=1:length(tspike)
        if newt_[i] == 1
            tspike[i] = t * dt
        end
    end
    # Reset membrane potential to resting potential after spike
    map!(x -> x >= v_thr ? EL : x, v, v)

    return v, Ibg, newt_, tspike
end
