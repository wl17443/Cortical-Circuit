using DelimitedFiles

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/logs/parameter_sweep/k_excexc/"
log_files = readdir(path)

reg_bump = r"Bump at neuron (?<mean>[+-]?([0-9]*[.])?[0-9]+), spreading across (?<sd>\d+) neurons."
reg_param = r"(?<param>\d+).txt"
reg_checkpoint = r"Checkpoint (?<checkpoint>\d+):"

# displacements = [[] for i in 1:length(log_files)]
# displacements = zeros(length(log_files), 4)
# param = 1

for file_nr=1:length(log_files)
    param = match(reg_param, log_files[file_nr])[:param]

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
                # displacements[file_nr, 4] = sum(abs.(diff(means)))
                map!(x -> x < 5 || x > 25 ? "E" : x, spreads, spreads)
                
                means = []
                spreads = []
                line_count = 1
            else
                line_count += 1
            end
        end
    end
end

results = open("results.csv", "w")
writedlm(results, displacements)
close(results)
