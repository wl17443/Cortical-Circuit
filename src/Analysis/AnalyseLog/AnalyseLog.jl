using DelimitedFiles

path = "C:/Users/Orion/Documents/Cortical-Circuit/logs/parameter_sweep/s_pvsst/"
log_files = readdir(path)

reg_bump = r"Bump at neuron (?<mean>[+-]?([0-9]*[.])?[0-9]+), spreading across (?<sd>\d+) neurons."
reg_param = r"(?<param>[+-]?([0-9]*[.])?[0-9]+).txt"
reg_checkpoint = r"Checkpoint (?<checkpoint>\d+):"

result = [[] for i in 1:length(log_files)]

for file_nr=1:length(log_files)
    param = match(reg_param, log_files[file_nr])[:param]
    push!(result[file_nr], param)
    open(path * log_files[file_nr]) do file
        line_count = 1
        means = []
        spreads = []
        for line in eachline(file)
            if line_count < 21
                m = match(reg_bump, line)
                if m != nothing
                    append!(means, parse(Float64, m[:mean]))
                    append!(spreads, parse(Float64, m[:sd]))
                end
                line_count += 1
            elseif line_count == 21
                m = match(reg_bump, line)
                if m != nothing
                    append!(means, parse(Float64, m[:mean]))
                    append!(spreads, parse(Float64, m[:sd]))
                end
                c = 0
                map(x -> x < 5 || x > 25 ? c += 1 : x, spreads)
                push!(result[file_nr], (sum(abs.(diff(means))), c))

                means = []
                spreads = []
                line_count = 1
            else
                line_count += 1
            end
        end
    end
end

results = open("s_pvsst.csv", "w")
writedlm(results, result)
close(results)
