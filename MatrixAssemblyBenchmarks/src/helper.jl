function get_folder_name(params::@NamedTuple{nruns::Int64, cells_per_dir::NTuple{N,Int64}, parts_per_dir::NTuple{N,Int64}}, root_name="") where {N}
    nruns, cells_per_dir, parts_per_dir = params
    cells_per_dir_str = sprint(show, cells_per_dir)
    cells_per_dir_str = replace(cells_per_dir_str, " " => "")
    parts_per_dir_str = sprint(show, parts_per_dir)
    parts_per_dir_str = replace(parts_per_dir_str, " " => "")
    nruns_str = sprint(show, nruns)
    mkpath(joinpath(root_name, join([cells_per_dir_str, parts_per_dir_str, nruns_str], "_")))
end

function get_folder_name(job_params, root_name="")
    nruns, cells_per_dir, parts_per_dir, _ = job_params
    cells_per_dir_str = sprint(show, cells_per_dir)
    cells_per_dir_str = replace(cells_per_dir_str, " " => "")
    parts_per_dir_str = sprint(show, parts_per_dir)
    parts_per_dir_str = replace(parts_per_dir_str, " " => "")
    nruns_str = sprint(show, nruns)
    mkpath(joinpath(root_name, join([cells_per_dir_str, parts_per_dir_str, nruns_str], "_")))
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

function create_job_file(node, core; project=Base.active_project(), dir_name="", cells_per_dirs=nothing, nrunss=nothing, methods=nothing)
    header_params = Dict(string(:node) => node, string(:core) => core)
    np = node * core
    project = abspath(project)
    parts_per_dir = node_core_partitions[(node, core)]
    root_name = "\"($node,$core)\""
    body_head_params = Dict(string(:np) => np, string(:project) => project, string(:parts_per_dir) => parts_per_dir, string(:root_name) => root_name)
    dir_name = mkpath(abspath(dir_name))
    file_name = joinpath(dir_name, "$(hash((header_params["node"], header_params["core"]))).sh")
    open(file_name, "w") do io
        render(io, template_header, header_params)
        if isnothing(cells_per_dirs)
            render(io, template_experiments_set, body_head_params)
        else
            render(io, template_body_head, body_head_params)
            if isnothing(methods)
                if eltype(cells_per_dirs) <: Tuple
                    for (cells_per_dir, nruns) in zip(cells_per_dirs, nrunss)
                        body_params = Dict(string(:cells_per_dir) => cells_per_dir, string(:nruns) => nruns)
                        render(io, template_experiements_body, body_params)
                    end
                else
                    body_params = Dict(string(:cells_per_dir) => cells_per_dirs, string(:nruns) => nrunss)
                    render(io, template_experiements_body, body_params)
                end
            else
                if eltype(cells_per_dirs) <: Tuple
                    for (cells_per_dir, nruns, method) in zip(cells_per_dirs, nrunss, methods)
                        body_params = Dict(string(:cells_per_dir) => cells_per_dir, string(:nruns) => nruns, string(:method) => "\"$method\"")
                        render(io, template_experiement_body, body_params)
                    end
                else
                    body_params = Dict(string(:cells_per_dir) => cells_per_dirs, string(:nruns) => nrunss)
                    render(io, template_experiements_body, body_params)
                end
            end
            write(io, template_body_tail)
        end
    end
    return file_name
end

function run_experiments_set(; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    files = Dict(create_job_file(node, core; project=project, dir_name=dir_name) => nexec for (node, core) in keys(node_core_partitions))
    cd(mkpath(abspath(dir_name)))
    while isempty(files)
        for file_name in keys(files)
            run(`$cmd $file_name`)
            files[file_name] -= 1
            if files[file_name] == 0
                delete!(files, file_name)
            end
        end
    end
end

function run_experiments(node, core, cells_per_dirs, nrunss; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    file_name = create_job_file(node, core; project=project, dir_name=dir_name, cells_per_dirs=cells_per_dirs, nrunss=nrunss)
    cd(mkpath(abspath(dir_name)))
    while nexec > 0
        run(`$cmd $file_name`)
        nexec -= 1
    end
end

function run_experiment(node, core, cells_per_dirs, nrunss, methods; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    file_name = create_job_file(node, core; project=project, dir_name=dir_name, cells_per_dirs=cells_per_dirs, nrunss=nrunss, methods=methods)
    cd(mkpath(abspath(dir_name)))
    while nexec > 0
        run(`$cmd $file_name`)
        nexec -= 1
    end
end
