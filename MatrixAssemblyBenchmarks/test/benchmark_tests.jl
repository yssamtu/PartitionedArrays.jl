import MatrixAssemblyBenchmarks as mb
dir_name = abspath("each_part")
nexec = 1
mb.run_experiments_sets(; dir_name=dir_name, nexec=nexec)
# nodes = Int[]
# cores = Int[]
# cells_per_dirss = NTuple{3, Int}[]
# nrunss = Int[]
# methods = String[]

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_auto_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 18)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_auto_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_no_compressed_snd_and_with_auto_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_int_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache")

# append!(nodes, 27)
# append!(cores, 1)
# push!(cells_per_dirss, (320, 320, 320))
# append!(nrunss, 20)
# push!(methods, "assemble_matrix_with_compressed_snd_and_with_auto_cache")

# append!(nodes, 18)
# append!(cores, 16)
# push!(cells_per_dirss, (80, 80, 80))
# append!(nrunss, 80)
# push!(methods, "petsc_coo")

# append!(nodes, 18)
# append!(cores, 16)
# push!(cells_per_dirss, (20, 20, 20))
# append!(nrunss, 320)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 12)
# push!(cells_per_dirss, (80, 80, 80))
# append!(nrunss, 80)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 12)
# push!(cells_per_dirss, (40, 40, 40))
# append!(nrunss, 160)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 12)
# push!(cells_per_dirss, (20, 20, 20))
# append!(nrunss, 320)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 16)
# push!(cells_per_dirss, (80, 80, 80))
# append!(nrunss, 80)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 16)
# push!(cells_per_dirss, (40, 40, 40))
# append!(nrunss, 160)
# push!(methods, "petsc_coo")

# append!(nodes, 27)
# append!(cores, 16)
# push!(cells_per_dirss, (20, 20, 20))
# append!(nrunss, 320)
# push!(methods, "petsc_coo")

# Base.exit_on_sigint(false)
# tasks = []
# @sync for (node, core) in zip(nodes, cores)
#     try
#         while length(tasks) >= 1
#             for i in length(tasks):-1:1
#                 if istaskdone(tasks[i])
#                     deleteat!(tasks, i)
#                 end
#             end
#             yield()
#         end

#         task = @async mb.run_experiments_set(node, core; dir_name=dir_name, nexec=nexec)
#         push!(tasks, task)
#     catch
#         break
#     end
# end

# Base.exit_on_sigint(false)
# tasks = []
# @sync for (node, core, cells_per_dir, nruns, method) in zip(nodes, cores, cells_per_dirss, nrunss, methods)
#     try
#         while length(tasks) >= 3
#             for i in length(tasks):-1:1
#                 if istaskdone(tasks[i])
#                     deleteat!(tasks, i)
#                 end
#             end
#             yield()
#         end

#         task = @async mb.run_experiment(node, core, cells_per_dir, nruns, method; dir_name=dir_name, nexec=nexec)
#         push!(tasks, task)
#     catch
#         break
#     end
# end
# mb.merge_dir("/home/ppp23002/thesis/strong_scaling/(18,12)", "/home/ppp23002/thesis/test/result/(18,12)")
# node=2
# core=1
# cells_per_dirs=(2, 1, 1)
# nrunss=1
# mb.run_experiments_sets(; dir_name=dir_name, nexec=1)
# mb.run_experiments_set(node, core; dir_name=dir_name)
# mb.run_experiments(node, core, cells_per_dirs, nrunss; dir_name=dir_name, nexec=1)
# mb.run_experiment(node, core, cells_per_dirs, nrunss, "psparse"; dir_name=dir_name)
# #(1, 1)
# parts_per_dir = (1, 1, 1)
# mb.experiments_set(parts_per_dir, "(1,1)")

# #(1, 2)
# parts_per_dir = (1, 1, 2)
# mb.experiments_set(parts_per_dir, "(1,2)")

# #(1, 4)
# parts_per_dir = (1, 2, 2)
# mb.experiments_set(parts_per_dir, "(1,4)")

# #(1, 8)
# parts_per_dir = (2, 2, 2)
# mb.experiments_set(parts_per_dir, "(1,8)")

# #(1, 12)
# parts_per_dir = (2, 2, 3)
# mb.experiments_set(parts_per_dir, "(1,12)")

# #(1, 16)
# parts_per_dir = (2, 2, 4)
# mb.experiments_set(parts_per_dir, "(1,16)")

# #(2, 1)
# parts_per_dir = (1, 1, 2)
# mb.experiments_set(parts_per_dir, "(2,1)")

# #(2, 2)
# parts_per_dir = (1, 1, 4)
# mb.experiments_set(parts_per_dir, "(2,2)")

# #(2, 4)
# parts_per_dir = (1, 2, 4)
# mb.experiments_set(parts_per_dir, "(2,4)")

# #(2, 8)
# parts_per_dir = (2, 2, 4)
# mb.experiments_set(parts_per_dir, "(2,8)")

# #(2, 12)
# parts_per_dir = (2, 2, 6)
# mb.experiments_set(parts_per_dir, "(2,12)")

# #(2, 16)
# parts_per_dir = (2, 2, 8)
# mb.experiments_set(parts_per_dir, "(2,16)")

# #(4, 1)
# parts_per_dir = (1, 2, 2)
# mb.experiments_set(parts_per_dir, "(4,1)")

# #(4, 2)
# parts_per_dir = (1, 2, 4)
# mb.experiments_set(parts_per_dir, "(4,2)")

# #(4, 4)
# parts_per_dir = (1, 4, 4)
# mb.experiments_set(parts_per_dir, "(4,4)")

# #(4, 8)
# parts_per_dir = (2, 4, 4)
# mb.experiments_set(parts_per_dir, "(4,8)")

# #(4, 12)
# parts_per_dir = (2, 4, 6)
# mb.experiments_set(parts_per_dir, "(4,12)")

# #(4, 16)
# parts_per_dir = (2, 4, 8)
# mb.experiments_set(parts_per_dir, "(4,16)")

# #(8, 1)
# parts_per_dir = (2, 2, 2)
# mb.experiments_set(parts_per_dir, "(8,1)")

# #(8, 2)
# parts_per_dir = (2, 2, 4)
# mb.experiments_set(parts_per_dir, "(8,2)")

# #(8, 4)
# parts_per_dir = (2, 4, 4)
# mb.experiments_set(parts_per_dir, "(8,4)")

# #(8, 8)
# parts_per_dir = (4, 4, 4)
# mb.experiments_set(parts_per_dir, "(8,8)")

# #(8, 12)
# parts_per_dir = (4, 4, 6)
# mb.experiments_set(parts_per_dir, "(8,12)")

# #(8, 16)
# parts_per_dir = (4, 4, 8)
# mb.experiments_set(parts_per_dir, "(8,16)")

# #(12, 1)
# parts_per_dir = (2, 2, 3)
# mb.experiments_set(parts_per_dir, "(12,1)")

# #(12, 2)
# parts_per_dir = (2, 2, 6)
# mb.experiments_set(parts_per_dir, "(12,2)")

# #(12, 4)
# parts_per_dir = (2, 4, 6)
# mb.experiments_set(parts_per_dir, "(12,4)")

# #(12, 8)
# parts_per_dir = (4, 4, 6)
# mb.experiments_set(parts_per_dir, "(12,8)")

# #(12, 12)
# parts_per_dir = (4, 4, 9)
# mb.experiments_set(parts_per_dir, "(12,12)")

# #(12, 16)
# parts_per_dir = (4, 4, 12)
# mb.experiments_set(parts_per_dir, "(12,16)")

# #(18, 1)
# parts_per_dir = (2, 3, 3)
# mb.experiments_set(parts_per_dir, "(18,1)")

# #(18, 2)
# parts_per_dir = (2, 3, 6)
# mb.experiments_set(parts_per_dir, "(18,2)")

# #(18, 4)
# parts_per_dir = (2, 6, 6)
# mb.experiments_set(parts_per_dir, "(18,4)")

# #(18, 8)
# parts_per_dir = (4, 6, 6)
# mb.experiments_set(parts_per_dir, "(18,8)")

# #(18, 12)
# parts_per_dir = (4, 6, 9)
# mb.experiments_set(parts_per_dir, "(18,12)")

# #(18, 16)
# parts_per_dir = (4, 6, 12)
# mb.experiments_set(parts_per_dir, "(18,16)")

# #(27, 1)
# parts_per_dir = (3, 3, 3)
# mb.experiments_set(parts_per_dir, "(27,1)")

# #(27, 2)
# parts_per_dir = (3, 3, 6)
# mb.experiments_set(parts_per_dir, "(27,2)")

# #(27, 4)
# parts_per_dir = (3, 6, 6)
# mb.experiments_set(parts_per_dir, "(27,4)")

# #(27, 8)
# parts_per_dir = (6, 6, 6)
# mb.experiments_set(parts_per_dir, "(27,8)")

# #(27, 12)
# parts_per_dir = (6, 6, 9)
# mb.experiments_set(parts_per_dir, "(27,12)")

# #(27, 16)
# parts_per_dir = (6, 6, 12)
# mb.experiments_set(parts_per_dir, "(27,16)")






# # N = 1
# cells_per_dir = (352, 352, 352)
# parts_per_dir = (1, 1, 1)
# nruns = 10

# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (88, 88, 88)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# # N = 2
# cells_per_dir = (352, 352, 352)
# parts_per_dir = (2, 1, 1)
# nruns = 10

# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (176, 88, 88)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# # N = 4
# cells_per_dir = (352, 352, 352)
# parts_per_dir = (2, 2, 1)
# nruns = 10

# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (176, 176, 88)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# # N = 8
# cells_per_dir = (176, 176, 176)
# parts_per_dir = (2, 2, 2)
# nruns = 10

# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (22, 22, 90112)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (22, 90112, 22)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (90112, 22, 22)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (22, 352, 5632)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (22, 5632, 352)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (352, 22, 5632)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (352, 5632, 22)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (5632, 22, 352)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (5632, 352, 22)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (352, 352, 352)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (1, 1, 8)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (1, 8, 1)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (8, 1, 1)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (1, 2, 4)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (1, 4, 2)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (2, 1, 4)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (2, 4, 1)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (4, 1, 2)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# parts_per_dir = (4, 2, 1)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# # N = 16
# cells_per_dir = (352, 352, 352)
# parts_per_dir = (4, 2, 2)
# nruns = 10

# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (352, 176, 176)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# N = 32
# cells_per_dir = (352, 352, 352)
# parts_per_dir = (4, 4, 2)
# nruns = 10

# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)

# cells_per_dir = (352, 352, 176)
# params = (; nruns, cells_per_dir, parts_per_dir)
# mb.experiments(params)
