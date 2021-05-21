#=============================================================================================================

Function to simulate two-compartmental pyramidal neurons. 

=============================================================================================================#

v_thr = -37.9e-3; v_rest = -68.1e-3

function simulatePyC1(t, dt, v_d, w_d, v_s, w_s, I_inj_d, t_, t_soma, st_SSTE, st_PVE, st_EE, st_E2E1, W_SSTE, W_PVE, W_EE, D)

    v_d += dv_d_dt(v_d, I_inj_d,  w_d, st_SSTE, W_SSTE, t_soma, t, D, dt) * dt
    w_d += dw_d_dt(w_d, v_d) * dt

    v_s += dv1_s_dt(v_s, v_d, w_s, st_PVE, st_EE, st_E2E1, W_PVE, W_EE, D, dt) * dt
    w_s += dw_s_dt(w_s, t_) * dt

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

    return v_d, w_d, v_s, w_s, newt_, t_soma
end

function simulatePyC2(t, dt, v_d, w_d, v_s, w_s, I_inj_s, t_, t_soma, st_SSTE, st_PVE, st_EE, st_E1E2, W_SSTE, W_PVE, W_EE, D)

    v_d += dv_d_dt(v_d, w_d, st_E1E2, st_SSTE, W_EE, W_SSTE, t_soma, t, D, dt) * dt
    w_d += dw_d_dt(w_d, v_d) .* dt

    v_s += dv2_s_dt(v_s, v_d, I_inj_s, w_s, st_PVE, st_EE, W_PVE, W_EE, D, dt) * dt
    w_s += dw_s_dt(w_s, t_) .* dt

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

    return v_d, w_d, v_s, w_s, newt_, t_soma
end
