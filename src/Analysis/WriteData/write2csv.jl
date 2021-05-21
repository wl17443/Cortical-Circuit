#=============================================================================================================

Functions used for writing membrane potential results to CSV. 

=============================================================================================================#

using CSV, DataFrames

function write2csv(arr::Array{Float64,2}, folder_name::String, datetime::String)
    csv_filename = datetime * "-voltage.csv"
    # params_filename = datetime * "-params.txt"

    csv_filepath = "C:/Users/Orion/Documents/University/Dissertation/Julia/data/$folder_name/$csv_filename"
    # params_filepath = "C:/Users/Orion/Documents/University/Dissertation/Julia/data/$folder_name/$params_filename"

    # CSV.write(params_filepath, con_params, header=false)
    CSV.write(csv_filepath, DataFrame(arr), header=false)
end
