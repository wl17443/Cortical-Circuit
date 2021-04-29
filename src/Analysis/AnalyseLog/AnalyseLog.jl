using DelimitedFiles

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/logs/noise_base_parameters/"
log_files = readdir(path)

reg_bump = r"Bump at neuron (?<mean>[+-]?([0-9]*[.])?[0-9]+), spreading across (?<sd>\d+) neurons."
reg_param = r"(?<som>\d+)-(?<den>\d+)-(?<inh>\d+).txt"
reg_checkpoint = r"Checkpoint (?<checkpoint>\d+):"

# displacements = [[] for i in 1:length(log_files)]
displacements = zeros(length(log_files), 4)
# param = 1

for file_nr=1:length(log_files)
    params = match(reg_param, log_files[file_nr])

    displacements[file_nr, 1] = parse(Int64, params[:som])
    displacements[file_nr, 2] = parse(Int64, params[:den])
    displacements[file_nr, 3] = parse(Int64, params[:inh])

    open(path * log_files[file_nr]) do file
        line_count = 1
        means = []
        for line in eachline(file)
            # if line_count > 8 && line_count < 28
            if line_count < 20
                m = match(reg_bump, line)
                if m != nothing
                    append!(means, parse(Float64, m[:mean]))
                elseif match(reg_checkpoint, line) == nothing
                    append!(means, 0)
                end
                line_count += 1
            elseif line_count == 20
                m = match(reg_bump, line)
                if m != nothing
                    append!(means, parse(Float64, m[:mean]))
                elseif match(reg_checkpoint, line) == nothing
                    append!(means, 0)
                end
                displacements[file_nr, 4] = sum(abs.(diff(means)))
                means = []
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
