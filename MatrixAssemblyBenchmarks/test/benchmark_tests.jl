import MatrixAssemblyBenchmarks as mb

cells_per_dir = (100, 50, 50)
parts_per_dir = (1, 1, 1)
nruns = 10

params = (; nruns, cells_per_dir, parts_per_dir)
mb.experiments(params)
