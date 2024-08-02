function get_folder_name(params::@NamedTuple{nruns::Int64, cells_per_dir::NTuple{N,Int64}, parts_per_dir::NTuple{N,Int64}}, root_name="") where {N}
    nruns, cells_per_dir, parts_per_dir = params
    cells_per_dir_str = sprint(show, cells_per_dir)
    cells_per_dir_str = replace(cells_per_dir_str, " " => "")
    parts_per_dir_str = sprint(show, parts_per_dir)
    parts_per_dir_str = replace(parts_per_dir_str, " " => "")
    nruns_str = sprint(show, nruns)
    mkpath(abspath(joinpath(root_name, join([cells_per_dir_str, parts_per_dir_str, nruns_str], "_"))))
end

function get_folder_name(job_params, root_name="")
    nruns, cells_per_dir, parts_per_dir, _ = job_params
    cells_per_dir_str = sprint(show, cells_per_dir)
    cells_per_dir_str = replace(cells_per_dir_str, " " => "")
    parts_per_dir_str = sprint(show, parts_per_dir)
    parts_per_dir_str = replace(parts_per_dir_str, " " => "")
    nruns_str = sprint(show, nruns)
    mkpath(abspath(joinpath(root_name, join([cells_per_dir_str, parts_per_dir_str, nruns_str], "_"))))
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

function create_job_file(node, core, experiment_type; project=Base.active_project(), dir_name="", cells_per_dirs=nothing, nrunss=nothing, methods=nothing)
    header_params = Dict(string(:node) => node, string(:core) => core)
    np = node * core
    project = abspath(project)
    parts_per_dir = node_core_partitions[(node, core)]
    root_name = "\"($node,$core)\""
    body_head_params = Dict(string(:np) => np, string(:project) => project, string(:parts_per_dir) => parts_per_dir, string(:root_name) => root_name, string(:experiment_type) => "\"$experiment_type\"")
    dir_name = mkpath(abspath(dir_name))
    if !isnothing(methods)
        file_name = joinpath(dir_name, "$(hash([node, core, cells_per_dirs, nrunss, methods])).sh")
    elseif !isnothing(cells_per_dirs)
        file_name = joinpath(dir_name, "$(hash([node, core, cells_per_dirs, nrunss])).sh")
    else
        file_name = joinpath(dir_name, "$(hash([node, core])).sh")
    end
    open(file_name, "w") do io
        render(io, template_header, header_params)
        if isnothing(cells_per_dirs)
            render(io, template_experiments_set, body_head_params)
        else
            render(io, template_body_head, body_head_params)
            if isnothing(methods)
                if eltype(cells_per_dirs) <: Tuple
                    for (cells_per_dir, nruns) in zip(cells_per_dirs, nrunss)
                        body_params = Dict(string(:cells_per_dir) => cells_per_dir, string(:nruns) => nruns, string(:experiment_type) => "\"$experiment_type\"")
                        render(io, template_experiements_body, body_params)
                    end
                else
                    body_params = Dict(string(:cells_per_dir) => cells_per_dirs, string(:nruns) => nrunss, string(:experiment_type) => "\"$experiment_type\"")
                    render(io, template_experiements_body, body_params)
                end
            else
                if eltype(cells_per_dirs) <: Tuple
                    for (cells_per_dir, nruns, method) in zip(cells_per_dirs, nrunss, methods)
                        body_params = Dict(string(:cells_per_dir) => cells_per_dir, string(:nruns) => nruns, string(:method) => "\"$method\"", string(:experiment_type) => "\"$experiment_type\"")
                        render(io, template_experiement_body, body_params)
                    end
                else
                    body_params = Dict(string(:cells_per_dir) => cells_per_dirs, string(:nruns) => nrunss, string(:method) => "\"$methods\"", string(:experiment_type) => "\"$experiment_type\"")
                    render(io, template_experiement_body, body_params)
                end
            end
            write(io, template_body_tail)
        end
    end
    return file_name
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
        new_book_path = joinpath(new_case, "$f.json")
        new_book = JSON.parsefile(new_book_path; dicttype=DataStructures.OrderedDict)
        result_book_path = joinpath(result_case, "$f.json")
        if isfile(result_book_path)
            result_book = JSON.parsefile(result_book_path; dicttype=DataStructures.OrderedDict)
        else
            result_book = new_book
        end
        for (time, new_data) in new_time_data
            if !haskey(result_summary, f)
                result_summary[f] = DataStructures.OrderedDict(time => new_data)
                summary_updated = true
                book_updated = true
            elseif !haskey(result_summary[f], time)
                result_summary[f][time] = new_data
                summary_updated = true
                book_updated = true
            elseif sum(new_data) < sum(result_summary[f][time])
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

function run_experiments_sets(experiment_type; dir_name="", project=Base.active_project(), nexec=1, parallel=20)
    tasks = []
    Base.exit_on_sigint(false)
    @sync for (node, core) in keys(node_core_partitions)
        try
            while length(tasks) >= parallel
                for i in length(tasks):-1:1
                    if istaskdone(tasks[i])
                        deleteat!(tasks, i)
                    end
                end
                yield()
            end
            task = @async run_experiments_set(node, core, experiment_type; dir_name=dir_name, project=project, nexec=nexec)
            push!(tasks, task)
        catch
            break
        end
    end
end

function run_experiments_set(node, core, experiment_type; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    original_dir = dirname(project)
    dir_name = joinpath(original_dir, basename(dir_name))
    file_name = create_job_file(node, core, experiment_type; project=project, dir_name=dir_name)
    cd(dir_name)
    result_dir = mkpath(joinpath(pwd(), "result"))
    data_path = joinpath(pwd(), "($node,$core)")
    Base.exit_on_sigint(false)
    while nexec > 0
        jobid = nothing
        try
            jobid = readchomp(`$cmd --parsable $file_name`)
            while true
                state = split(readchomp(`sacct --jobs=$jobid --noheader --format=state --parsable2`), "\n")[1]
                if state == "COMPLETED"
                    rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
                    break
                elseif state == "FAILED" || occursin("CANCELLED", state)
                    break
                end
            end
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec -= 1
        catch
            rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
            try
                run(pipeline(`scancel $jobid`; stdout=devnull, stderr=devnull))
            catch
            end
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec = 0
        end
    end
end

function run_experiments(node, core, cells_per_dirs, nrunss, experiment_type; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    original_dir = dirname(project)
    dir_name = joinpath(original_dir, basename(dir_name))
    file_name = create_job_file(node, core, experiment_type; project=project, dir_name=dir_name, cells_per_dirs=cells_per_dirs, nrunss=nrunss)
    cd(dir_name)
    result_dir = mkpath(joinpath(pwd(), "result"))
    data_path = joinpath(pwd(), "($node,$core)")
    Base.exit_on_sigint(false)
    while nexec > 0
        jobid = nothing
        try
            jobid = readchomp(`$cmd --parsable $file_name`)
            while true
                state = split(readchomp(`sacct --jobs=$jobid --noheader --format=state --parsable2`), "\n")[1]
                if state == "COMPLETED"
                    rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
                    break
                elseif state == "FAILED" || occursin("CANCELLED", state)
                    break
                end
            end
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec -= 1
        catch
            rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
            try
                run(pipeline(`scancel $jobid`; stdout=devnull, stderr=devnull))
            catch
            end
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec = 0
        end
    end
end

function run_experiment(node, core, cells_per_dirs, nrunss, methods, experiment_type; dir_name="", project=Base.active_project(), nexec=1)
    cmd = :sbatch
    original_dir = dirname(project)
    dir_name = joinpath(original_dir, basename(dir_name))
    file_name = create_job_file(node, core, experiment_type; project=project, dir_name=dir_name, cells_per_dirs=cells_per_dirs, nrunss=nrunss, methods=methods)
    cd(dir_name)
    result_dir = mkpath(joinpath(pwd(), "result"))
    data_path = joinpath(pwd(), "($node,$core)")
    Base.exit_on_sigint(false)
    while nexec > 0
        jobid = nothing
        try
            jobid = readchomp(`$cmd --parsable $file_name`)
            while true
                state = split(readchomp(`sacct --jobs=$jobid --noheader --format=state --parsable2`), "\n")[1]
                if state == "COMPLETED"
                    rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
                    break
                elseif state == "FAILED" || occursin("CANCELLED", state)
                    break
                end
            end
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec -= 1
        catch
            rm(joinpath(pwd(), "slurm-$jobid.out"); force=true)
            try
                run(pipeline(`scancel $jobid`; stdout=devnull, stderr=devnull))
            catch
            end
            result_path = joinpath(result_dir, basename(data_path))
            merge_dir(result_path, data_path)
            nexec = 0
        end
    end
end
