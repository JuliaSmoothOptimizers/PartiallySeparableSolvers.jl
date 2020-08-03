using CSV

repo_save = "src/comparaison/file_dataframe/"

function save_dataframe(stats)
    keys_stats = keys(stats)
    for i in keys_stats
        file = repo_save * string(i) * ".csv"
        CSV.write(file, stats[i])
    end
end

function get_dataframe()
    dic_dataframe = Dict{Symbol,DataFrame}()
    for i in keys_stats
        string_name = string(i)
        file = repo_save * string_name * ".csv"
        symbol_name = Symbol(string_name)
        df =  DataFrame(CSV.File( file ))
        dic_dataframe[symbol_name] = df
    end
end
