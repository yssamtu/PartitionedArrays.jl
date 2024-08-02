function get_execution_time(path)
    function get_maximum(v)
        np = length(v)
        nruns = length(v[1])
        while np > 1
            left = 1
            right = np
            while left < right
                v[left] = max.(v[left], v[right])
                left += 1
                right -= 1
            end
            np = (np + 1) >> 1
        end
        v[1]
    end
    json_dict = JSON.parsefile(path)
    buildmat = Vector{Vector{Float64}}(json_dict["buildmat"])
    rebuildmat = Vector{Vector{Float64}}(json_dict["rebuildmat"])
    build_time = minimum(get_maximum(buildmat))
    rebuild_time = minimum(get_maximum(rebuildmat))
    (; build_time, rebuild_time)
end

function get_execution_time(path, buildmat, rebuildmat)
    function get_maximum(v)
        np = length(v)
        nruns = length(v[1])
        while np > 1
            left = 1
            right = np
            while left < right
                v[left] = max.(v[left], v[right])
                left += 1
                right -= 1
            end
            np = (np + 1) >> 1
        end
        v[1]
    end
    build_time = minimum(get_maximum(buildmat))
    rebuild_time = minimum(get_maximum(rebuildmat))
    (; build_time, rebuild_time)
end

function get_all_execution_times(params)
    folder_name = get_folder_name(params)
    execution_times = DataStructures.OrderedDict{String,@NamedTuple{build_time::Float64, rebuild_time::Float64}}()
    for method in methods
        path = get_path(method, folder_name)
        execution_times[method] = get_execution_time(path)
    end
    execution_times
end
