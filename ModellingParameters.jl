module ModellingParameters

include("Units.jl")

using .Units 

## Model Constants 
nr_pyc = 1
nr_chc = 1
nr_sst = 1
nr_pv = 1

t = 100*ms
dt = 1*ms

steps = Int(t/dt)

S = 0.1

# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 