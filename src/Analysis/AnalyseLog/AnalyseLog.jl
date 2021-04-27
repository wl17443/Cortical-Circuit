using DelimitedFiles

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/logs/noise_base_parameters/"
log_files = readdir(path)

reg_bump = r"Bump at neuron (?<mean>[+-]?([0-9]*[.])?[0-9]+), spreading across (?<sd>\d+) neurons."
reg_param = r"(?<param>[+-]?([0-9]*[.])?[0-9]+e-[0-9])-(?<noise>\d+).txt"
reg_checkpoint = r"Checkpoint (?<checkpoint>\d+):"

displacements = [[] for i in 1:11]
# param = 1

for file_nr=1:length(log_files)
    # if file_nr % 10 == 1 && file_nr > 10 && file_nr < 11
    #     global param += 1
    # end
    open(path * log_files[file_nr]) do file
        line_count = 1
        means = []
        for line in eachline(file)
            if line_count > 8 && line_count < 28
                m = match(reg_bump, line)
                if m != nothing
                    append!(means, parse(Float64, m[:mean]))
                elseif match(reg_checkpoint, line) == nothing
                    append!(means, 0)
                end
                line_count += 1
            elseif line_count == 28
                m = match(reg_bump, line)
                if m != nothing
                    append!(means, parse(Float64, m[:mean]))
                elseif match(reg_checkpoint, line) == nothing
                    append!(means, 0)
                end
                append!(displacements[file_nr], sum(abs.(diff(means))))
                # println(log_files[file_nr])
                # println(sum(abs.(diff(means))))
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
