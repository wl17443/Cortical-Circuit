using DelimitedFiles

file_path = "C:/Users/Orion/Documents/University/Dissertation/Julia/logs/07-04-2021-18-39.txt"

disfunctional_idx = []
functional_idx    = []

open(file_path, "r") do file
    params_counter = 1
    status_counter = 1
    experiment = 1

    params = Dict()
    status = []
    for line in eachline(file)
        if params_counter <= 8
            pair = split(line, " ")
            params[pair[1]] = parse(Float64, pair[2])
            params_counter += 1
        else if status_counter <= 20
             push!(status, line)
             status_counter += 1
        else
            ## Analyse current parameter mode

            ## IO interactions
            # println("Analysing Experiment $experiment.")
            # println("Show heatmap? [y/n]")
            # response = readline()
            # if response == "y"
            # end


            ## Reset counters and variables
            params_counter = 1
            status_counter = 1
            params = Dict()
            status = []
            experiment += 1
        end
    end
end
