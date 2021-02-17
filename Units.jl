module Units 

mV = 0.001
ms = 0.001
pA = 1e-9
pF = 1e-9
nS = 1e-6
nA = 1e-6
s = 1

# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end