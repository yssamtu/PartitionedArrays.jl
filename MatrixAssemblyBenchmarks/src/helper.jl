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
                    body_params = Dict(string(:cells_per_dir) => cells_per_dirs, string(:nruns) => nrunss, string(:method) => "\"$methods\"")
                    render(io, template_experiement_body, body_params)
                end
            end
            write(io, template_body_tail)
        end
    end
    return file_name
end

function run_experiments_sets(; dir_name="", project=Base.active_project(), nexec=1)
    original_dir = pwd()
    dir_name = abspath(dir_name)
    @sync for (node, core) in keys(node_core_partitions)
        @async run_experiments_set(node, core; dir_name=dir_name, project=project, nexec=nexec)
    end
    cd(original_dir)
end

function run_experiments_set(node, core; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    file_name = create_job_file(node, core; project=project, dir_name=dir_name)
    original_dir = pwd()
    cd(abspath(dir_name))
    result_dir = mkpath(joinpath(pwd(), "result"))
    data_path = joinpath(pwd(), "($node,$core)")
    Base.exit_on_sigint(false)
    while nexec > 0
        jobid = nothing
        try
            jobid = readchomp(`$cmd --parsable $file_name`)
            while true
                state = split(readchomp(`sacct --jobs=$jobid --noheader --format=state --parsable2`), "\n")[1]
                if state == "COMPLETED" || state == "FAILED"
                    break
                end
            end
            rm(joinpath(pwd(), "slurm-$jobid.out"))
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec -= 1
        catch
            rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
            try
                run(pipeline(`preserve -c $jobid`; stdout=devnull, stderr=devnull))
            catch
            end
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec = 0
        end
    end
    Base.exit_on_sigint(true)
    cd(original_dir)
end

function run_experiments(node, core, cells_per_dirs, nrunss; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    file_name = create_job_file(node, core; project=project, dir_name=dir_name, cells_per_dirs=cells_per_dirs, nrunss=nrunss)
    original_dir = pwd()
    cd(abspath(dir_name))
    result_dir = mkpath(joinpath(pwd(), "result"))
    data_path = joinpath(pwd(), "($node,$core)")
    Base.exit_on_sigint(false)
    while nexec > 0
        jobid = nothing
        try
            jobid = readchomp(`$cmd --parsable $file_name`)
            while true
                state = split(readchomp(`sacct --jobs=$jobid --noheader --format=state --parsable2`), "\n")[1]
                if state == "COMPLETED" || state == "FAILED"
                    break
                end
            end
            rm(joinpath(pwd(), "slurm-$jobid.out"))
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec -= 1
        catch
            rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
            run(pipeline(`preserve -c $jobid`, devnull))
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec = 0
        end
    end
    Base.exit_on_sigint(true)
    cd(original_dir)
end

function run_experiment(node, core, cells_per_dirs, nrunss, methods; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    file_name = create_job_file(node, core; project=project, dir_name=dir_name, cells_per_dirs=cells_per_dirs, nrunss=nrunss, methods=methods)
    original_dir = pwd()
    cd(abspath(dir_name))
    result_dir = mkpath(joinpath(pwd(), "result"))
    data_path = joinpath(pwd(), "($node,$core)")
    Base.exit_on_sigint(false)
    while nexec > 0
        jobid = nothing
        try
            jobid = readchomp(`$cmd --parsable $file_name`)
            while true
                state = split(readchomp(`sacct --jobs=$jobid --noheader --format=state --parsable2`), "\n")[1]
                if state == "COMPLETED" || state == "FAILED"
                    break
                end
            end
            rm(joinpath(pwd(), "slurm-$jobid.out"))
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec -= 1
        catch
            rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
            run(pipeline(`preserve -c $jobid`, devnull))
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec = 0
        end
    end
    Base.exit_on_sigint(true)
    cd(original_dir)
end

function merge_file(result_case, new_case)
    new_summary_path = joinpath(new_case, "summary.json")
    if !isfile(new_summary_path)
        return
    end
    new_summary = JSON.parsefile(new_summary_path; dicttype=DataStructures.OrderedDict)
    result_summary_path = joinpath(result_case, "summary.json")
    result_summary = JSON.parsefile(result_summary_path; dicttype=DataStructures.OrderedDict)
    summary_updated = false
    book_updated = false
    for (f, new_time_data) in new_summary
        result_book_path = joinpath(result_case, "$f.json")
        result_book = JSON.parsefile(result_book_path; dicttype=DataStructures.OrderedDict)
        new_book_path = joinpath(new_case, "$f.json")
        new_book = JSON.parsefile(new_book_path; dicttype=DataStructures.OrderedDict)
        for (time, new_data) in new_time_data
            if new_data < result_summary[f][time]
                result_summary[f][time] = new_data
                k = "$(time[1:end-5])mat"
                result_book[k] = new_book[k]
                summary_updated = true
                book_updated = true
            end
        end
        if book_updated
            open(result_book_path, "w") do f
                JSON.print(f, result_book, 2)
            end
            book_updated = false
        end
        rm(new_book_path)
    end
    if summary_updated
        open(result_summary_path, "w") do f
            JSON.print(f, result_summary, 2)
        end
    end
    rm(new_summary_path)
end

function merge_dir(result_dir, new_dir)
    if !isdir(new_dir)
        return
    end
    new_cases = readdir(new_dir, join=true)
    for new_case in new_cases
        result_case = joinpath(result_dir, basename(new_case))
        if isdir(result_case)
            merge_file(result_case, new_case)
            rm(new_case; recursive=true)
        else
            new_summary_path = joinpath(new_case, "summary.json")
            if isfile(new_summary_path)
                mkpath(result_dir)
                mv(new_case, result_case)
            else
                rm(new_case; recursive=true)
            end
        end
    end
    rm(new_dir)
end
