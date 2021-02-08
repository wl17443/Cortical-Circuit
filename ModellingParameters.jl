module ModellingParameters

include("Units.jl")

using .Units 

## Model Constants 
nr_pyc = 5
nr_chc = 1
nr_sst = 5
nr_pv = 5

t = 100*s
dt = 5*ms

# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end 