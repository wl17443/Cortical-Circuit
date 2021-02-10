module AxonInitialSegment

include("Units.jl")

using .Units 

EL = -70*mV; t_a = 14*ms ## TODO - Subject to change 
g_s = 1000*pA ## TODO - Subject to change 

dv_a_dt(v_a, t) = -(v_a .- EL) ./ t_a + (g_s * I_s + I_chc) ./ C_a

end 