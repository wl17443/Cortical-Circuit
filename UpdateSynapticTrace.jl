module UpdateSynapticTrace 

include("Units.jl")
include("ModellingParameters.jl")

using .Units
using .ModellingParameters

function update_st(t, t_, st)
    if t_ == t*dt
        st += S
    else 
        st = st + ds_dt(st)*dt
    end 

    return st
end 

t_s = 5*ms

function ds_dt(st)
    return - st / t_s
end 

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 