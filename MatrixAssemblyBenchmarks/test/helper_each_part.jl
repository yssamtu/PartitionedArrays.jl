function get_folder_name(params::@NamedTuple{nruns::Int64, cells_per_dir::NTuple{N, Int64}, parts_per_dir::NTuple{N, Int64}}) where {N}
    nruns, cells_per_dir, parts_per_dir = params
    cells_per_dir_str = sprint(show, cells_per_dir)
    cells_per_dir_str = replace(cells_per_dir_str, " " => "")
    parts_per_dir_str = sprint(show, parts_per_dir)
    parts_per_dir_str = replace(parts_per_dir_str, " " => "")
    nruns_str = sprint(show, nruns)
    mkpath(join([cells_per_dir_str, parts_per_dir_str, nruns_str], "_"))
end

function get_folder_name(job_params)
    nruns, cells_per_dir, parts_per_dir, _ = job_params
    cells_per_dir_str = sprint(show, cells_per_dir)
    cells_per_dir_str = replace(cells_per_dir_str, " " => "")
    parts_per_dir_str = sprint(show, parts_per_dir)
    parts_per_dir_str = replace(parts_per_dir_str, " " => "")
    nruns_str = sprint(show, nruns)
    mkpath(join([cells_per_dir_str, parts_per_dir_str, nruns_str], "_"))
end

function get_path(job_params, folder_name=get_folder_name(job_params))
    file_name = job_params.method
    extension = ".json"
    joinpath(folder_name, file_name * extension)
end

function get_path(method::String, folder_name)
    file_name = method
    extension = ".json"
    joinpath(folder_name, file_name * extension)
end

function get_params(path)
    function parse_tuple(type, str)
        vec_str = split(str[2:end-1], ",")
        vec_type = parse.(type, vec_str)
        NTuple{length(vec_type),eltype(vec_type)}(vec_type)
    end
    abs_path = abspath(path)
    dir_name = basename(dirname(abs_path))
    base_name = basename(abs_path)
    cells_per_dir_str, parts_per_dir_str, nruns_str = split(dir_name, "_")
    nruns = parse(Int, nruns_str)
    cells_per_dir = parse_tuple(Int, cells_per_dir_str)
    parts_per_dir = parse_tuple(Int, parts_per_dir_str)
    method = splitext(base_name)[1]
    (; nruns, cells_per_dir, parts_per_dir, method)
end

function get_execution_time(path)
    function get_maximum(v)
        np = length(v)
        nruns = length(v[1])
        while np > 1
            left = 1
            right = np
            while left < right
                for i in 1:nruns
                    if v[left][i][1]+v[left][i][2]+v[left][i][3] < v[right][i][1]+v[right][i][2]+v[right][i][3]
                        v[left] = v[right]
                    end
                end
                # v[left] = max.(v[left], v[right])
                left += 1
                right -= 1
            end
            np = (np + 1) >> 1
        end
        v[1]
    end
    json_dict = JSON.parsefile(path)
    buildmat = Vector{Vector{Pair{Float64,Tuple{Float64,Float64}}}}(json_dict["buildmat"])
    rebuildmat = Vector{Vector{Pair{Float64,Tuple{Float64,Float64}}}}(json_dict["rebuildmat"])
    max_build_time = get_maximum(buildmat)
    max_rebuild_time = get_maximum(rebuildmat)
    build_time = max_build_time[1]
    rebuild_time = max_rebuild_time[1]
    for i in 2:length(buildmat[1])
        if max_build_time[i][1]+max_build_time[i][2]+max_build_time[i][3] < build_time[1]+build_time[2]+build_time[3]
            build_time = max_build_time[i]
        end
        if max_rebuild_time[i][1]+max_rebuild_time[i][2]+max_rebuild_time[i][3] < rebuild_time[1]+rebuild_time[2]+rebuild_time[3]
            rebuild_time = max_rebuild_time[i]
        end
    end
    # build_time = minimum(get_maximum(buildmat))
    # rebuild_time = minimum(get_maximum(rebuildmat))
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
                for i in 1:nruns
                    if v[left][i][1]+v[left][i][2]+v[left][i][3] < v[right][i][1]+v[right][i][2]+v[right][i][3]
                        v[left] = v[right]
                    end
                end
                # v[left] = max.(v[left], v[right])
                left += 1
                right -= 1
            end
            np = (np + 1) >> 1
        end
        v[1]
    end
    max_build_time = get_maximum(buildmat)
    max_rebuild_time = get_maximum(rebuildmat)
    build_time = max_build_time[1]
    rebuild_time = max_rebuild_time[1]
    for i in 2:length(buildmat[1])
        if max_build_time[i][1]+max_build_time[i][2]+max_build_time[i][3] < build_time[1]+build_time[2]+build_time[3]
            build_time = max_build_time[i]
        end
        if max_rebuild_time[i][1]+max_rebuild_time[i][2]+max_rebuild_time[i][3] < rebuild_time[1]+rebuild_time[2]+rebuild_time[3]
            rebuild_time = max_rebuild_time[i]
        end
    end
    # build_time = minimum(get_maximum(buildmat))
    # rebuild_time = minimum(get_maximum(rebuildmat))
    (; build_time, rebuild_time)
end

function get_all_execution_times(params)
    folder_name = get_folder_name(params)
    execution_times = DataStructures.OrderedDict{String, @NamedTuple{build_time::Float64, rebuild_time::Float64}}()
    for method in methods
        path = get_path(method, folder_name)
        execution_times[method] = get_execution_time(path)
    end
    execution_times
end
