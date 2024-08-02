import MatrixAssemblyBenchmarks as mb
dir_name = abspath("test1")
experiment_type = "strong" # "strong", "weak", "each_part"
parallel = 1
nexec = 1
# mb.run_experiments_sets(experiment_type; dir_name=dir_name, nexec=nexec, parallel=parallel)
nodes = Int[]
cores = Int[]
cells_per_dirs = NTuple{3, Int}[]
nrunss = Int[]
methods = String[]

append!(nodes, 18)
append!(cores, 1)
push!(cells_per_dirs, (320, 320, 320))
append!(nrunss, 20)
push!(methods, "assemble_matrix_no_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_auto_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_auto_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_auto_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirs, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_auto_cache")

# append!(nodes, 18)
# append!(cores, 16)
# push!(cells_per_dirs, (80, 80, 80))
# append!(nrunss, 80)
# push!(methods, "petsc_coo")

# append!(nodes, 18)
# append!(cores, 16)
# push!(cells_per_dirs, (20, 20, 20))
# append!(nrunss, 320)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 12)
# push!(cells_per_dirs, (80, 80, 80))
# append!(nrunss, 80)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 12)
# push!(cells_per_dirs, (40, 40, 40))
# append!(nrunss, 160)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 12)
# push!(cells_per_dirs, (20, 20, 20))
# append!(nrunss, 320)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 16)
# push!(cells_per_dirs, (80, 80, 80))
# append!(nrunss, 80)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 16)
# push!(cells_per_dirs, (40, 40, 40))
# append!(nrunss, 160)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 16)
# push!(cells_per_dirs, (20, 20, 20))
# append!(nrunss, 320)
# push!(methods, "petsc_coo")

Base.exit_on_sigint(false)
tasks = []
# @sync for (node, core) in zip(nodes, cores)
@sync for (node, core, cells_per_dir, nruns, method) in zip(nodes, cores, cells_per_dirs, nrunss, methods)
    try
        while length(tasks) >= parallel
            for i in length(tasks):-1:1
                if istaskdone(tasks[i])
                    deleteat!(tasks, i)
                end
            end
            yield()
        end

        # task = @async mb.run_experiments_set(node, core, experiment_type; dir_name=dir_name, nexec=nexec)
        task = @async mb.run_experiment(node, core, cells_per_dir, nruns, method, experiment_type; dir_name=dir_name, nexec=nexec)
        push!(tasks, task)
    catch
        break
    end
end

# mb.merge_dir("/home/ppp23002/thesis/strong_scaling/(18,12)", "/home/ppp23002/thesis/test/result/(18,12)")
