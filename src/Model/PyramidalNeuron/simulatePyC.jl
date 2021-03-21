v_thr = -37.9e-3; v_rest = -68.1e-3

function simulatePyC(t, dt, v_d, w_d, v_s, w_s, I_inj_d, I_inj_s, I_dbg, I_sbg, bgnoise_lvl, t_, t_soma, st_SSTE, st_PVE, st_EE, W_SSTE, W_PVE, W_EE)

    v_d += dv_d_dt(v_d, I_inj_d, I_dbg, w_d, st_SSTE, W_SSTE, t_soma, t) .* dt
    w_d += dw_d_dt(w_d, v_d) .* dt

    v_s += dv_s_dt(v_s, v_d, I_inj_s, I_sbg, w_s, st_PVE, st_EE, W_PVE, W_EE) .* dt
    w_s += dw_s_dt(w_s, t_) .* dt

    I_dbg += dI_dbg_dt(I_dbg, bgnoise_lvl, dt) .* dt
    I_sbg += dI_sbg_dt(I_sbg, bgnoise_lvl, dt) .* dt

    newt_ = map(x -> x >= v_thr ? 1 : 0, v_s)
    for i=1:length(t_soma)
        if newt_[i] == 1
            t_soma[i] = t*dt
        end
    end
    map!(x -> x >= v_thr ? v_rest : x, v_s, v_s)

    return v_d, w_d, v_s, w_s, I_dbg, I_sbg, newt_, t_soma
end
