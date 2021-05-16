using NamedArrays
using DelimitedFiles

path = "C:/Users/Orion/Documents/University/Dissertation/Julia/logs/big_noise/s_excexc/"
log_files = readdir(path)

reg_bump = r"Bump at neuron (?<mean>[+-]?([0-9]*[.])?[0-9]+), spreading across (?<sd>\d+) neurons."
reg_param = r"(?<param>[+-]?([0-9]*[.])?[0-9]+)-(?<noise>\d+).txt"
reg_checkpoint = r"Checkpoint (?<checkpoint>\d+):"

results = NamedArray(zeros(6, 16))
setnames!(results, ["0", "1", "10", "100", "1000", "10000"], 1)
pars = collect(575.0:5.0:650.0)
setnames!(results, string.(pars), 2)

errors = NamedArray(zeros(6, 16))
setnames!(errors, ["0", "1", "10", "100", "1000", "10000"], 1)
pars = collect(575.0:5.0:650.0)
setnames!(errors, string.(pars), 2)

for file_nr=1:length(log_files)
    params = match(reg_param, log_files[file_nr])
    param = String(params[:param])
    noise = String(params[:noise])
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
                errors[noise, param] += c
                results[noise, param] += sum(abs.(diff(means)))

                means = []
                spreads = []
                line_count = 1
            else
                line_count += 1
            end
        end
    end
end

results ./= 5

result = open("big_noise_s_excexc.csv", "w")
writedlm(result, results)
close(result)

errors ./= 5

error = open("big_noise_s_excexc_error.csv", "w")
writedlm(error, errors)
close(error)
