#=============================================================================================================

Function to simulate two-compartmental pyramidal neurons

=============================================================================================================#

v_thr = -37.9e-3; v_rest = -68.1e-3

function simulatePyC(t, dt, v_d, w_d, v_s, w_s, noise_lvl, I_inj_d, I_inj_s, I_dbg, I_sbg, t_, t_soma, st_SSTE, st_PVE, st_EE, W_SSTE, W_PVE, W_EE)

    v_d += dv_d_dt(v_d, I_inj_d, I_dbg, w_d, st_SSTE, W_SSTE, t_soma, t) .* dt
    w_d += dw_d_dt(w_d, v_d) .* dt

    v_s += dv_s_dt(v_s, v_d, I_inj_s, I_sbg, w_s, st_PVE, st_EE, W_PVE, W_EE) .* dt
    w_s += dw_s_dt(w_s, t_) .* dt

    I_dbg += dI_dbg_dt(I_dbg, dt, noise_lvl) .* dt
    I_sbg += dI_sbg_dt(I_sbg, dt, noise_lvl) .* dt

    newt_ = map(x -> x >= v_thr ? 1 : 0, v_s)
    for i=1:length(t_soma)
        if t * dt - t_soma[i] <= 3e-3
            v_s[i] = v_rest
        end
        if newt_[i] == 1
            t_soma[i] = t * dt
            v_s[i] = v_rest
        end
    end

    return v_d, w_d, v_s, w_s, I_dbg, I_sbg, newt_, t_soma
end
