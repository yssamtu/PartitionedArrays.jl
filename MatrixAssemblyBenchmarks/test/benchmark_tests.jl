import MatrixAssemblyBenchmarks as mb
import PartitionedArrays as pa

cells_per_dir = (2, 2, 1)
parts_per_dir = (2, 2, 1)
nruns = 10

params = (; nruns, cells_per_dir, parts_per_dir, method="psparse")
out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

# params = (; nruns, cells_per_dir, parts_per_dir, method="petsc_setvalues")
# out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

# params = (; nruns, cells_per_dir, parts_per_dir, method="petsc_coo")
# out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_no_compressed_snd_and_with_int_vector_cache")
out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache")
out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_no_compressed_snd_and_with_auto_cache")
out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_with_compressed_snd_and_with_int_vector_cache")
out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache")
out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_with_compressed_snd_and_with_auto_cache")
out = pa.with_mpi(distribute -> mb.benchmark_psparse(distribute, params))


cells_per_dir = (500, 100, 100)
parts_per_dir = (2, 2, 1)
nruns = 10

params = (; nruns, cells_per_dir, parts_per_dir, method="psparse")
mb.experiment(params)

# params = (; nruns, cells_per_dir, parts_per_dir, method="petsc_coo")
# mb.experiment(params)

# params = (; nruns, cells_per_dir, parts_per_dir, method="petsc_setvalues")
# mb.experiment(params)

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_no_compressed_snd_and_with_int_vector_cache")
mb.experiment(params)

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache")
mb.experiment(params)

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_no_compressed_snd_and_with_auto_cache")
mb.experiment(params)

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_with_compressed_snd_and_with_int_vector_cache")
mb.experiment(params)

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache")
mb.experiment(params)

params = (; nruns, cells_per_dir, parts_per_dir, method="assemble_matrix_with_compressed_snd_and_with_auto_cache")
mb.experiment(params)