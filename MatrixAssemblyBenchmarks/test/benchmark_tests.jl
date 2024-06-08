import MatrixAssemblyBenchmarks as mb

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


cells_per_dir = (352, 352, 352)
parts_per_dir = (4, 2, 2)
nruns = 10

params = (; nruns, cells_per_dir, parts_per_dir)
mb.experiments(params)

cells_per_dir = (352, 176, 176)
params = (; nruns, cells_per_dir, parts_per_dir)
mb.experiments(params)
