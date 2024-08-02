methods = [
    "psparse",
    "petsc_setvalues",
    "petsc_coo",
    "assemble_matrix_no_compressed_snd_and_with_int_vector_cache",
    "assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache",
    "assemble_matrix_no_compressed_snd_and_with_auto_cache",
    "assemble_matrix_with_compressed_snd_and_with_int_vector_cache",
    "assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache",
    "assemble_matrix_with_compressed_snd_and_with_auto_cache"
]

if ENV == "das5"
    template_header = raw"""
    #!/bin/bash
    #SBATCH --time=00:15:00
    #SBATCH -N {{node}}
    #SBATCH --ntasks-per-node={{core}}

    module unload openmpi
    module unload intel/mkl
    module unload julia
    module load openmpi/gcc/64/4.0.2
    module load intel/mkl/64/11.2/2015.5.223
    module load julia/1.10.3

    """
else
    template_header = raw"""
    #!/bin/bash
    #SBATCH --time=06:00:00
    #SBATCH -N {{node}}
    #SBATCH --ntasks-per-node={{core}}

    module purge
    module load 2023
    module load imkl/2023.1.0
    module load OpenMPI/4.1.5-GCC-12.3.0
    module load juliaup/1.14.5-GCCcore-12.3.0
    juliaup add 1.10.4~x64
    juliaup default 1.10.4~x64

    """
end

template_experiments_set = raw"""
mpiexec -np {{np}} julia -O3 --check-bound=no --project={{{project}}} -e '
    import MatrixAssemblyBenchmarks as mb
    parts_per_dir = {{parts_per_dir}}
    experiment_type = {{{experiment_type}}}
    mb.experiments_set(parts_per_dir, experiment_type; root_name={{{root_name}}})
'
"""

template_body_head = raw"""
mpiexec -np {{np}} julia -O3 --check-bound=no --project={{{project}}} -e '
    import PartitionedArrays as pa
    import MatrixAssemblyBenchmarks as mb
    pa.with_mpi() do distribute
        root_name = {{{root_name}}}
        parts_per_dir = {{parts_per_dir}}
"""

template_experiements_body = raw"""
        cells_per_dir = {{cells_per_dir}}
        nruns = {{nruns}}
        params = (; nruns, cells_per_dir, parts_per_dir)
        experiment_type = {{{experiment_type}}}

        mb.experiments(params, experiment_type; root_name=root_name, distribute=distribute)
"""

template_experiement_body = raw"""
        cells_per_dir = {{cells_per_dir}}
        nruns = {{nruns}}
        method = {{{method}}}
        params = (; nruns, cells_per_dir, parts_per_dir, method)
        experiment_type = {{{experiment_type}}}
        
        mb.experiment(params, experiment_type; root_name=root_name, distribute=distribute)
"""

template_body_tail = raw"""
    end
'
"""

node_config = DataStructures.OrderedDict(1 => (1, 1, 1), 2 => (1, 1, 2), 4 => (1, 2, 2), 8 => (2, 2, 2), 12 => (2, 2, 3), 18 => (2, 3, 3), 27 => (3, 3, 3))
core_config = DataStructures.OrderedDict(1 => (1, 1, 1), 2 => (1, 1, 2), 4 => (1, 2, 2), 8 => (2, 2, 2), 12 => (2, 2, 3), 16 => (2, 2, 4))
# node_config = DataStructures.OrderedDict(27 => (3, 3, 3), 18 => (2, 3, 3), 12 => (2, 2, 3), 8 => (2, 2, 2), 4 => (1, 2, 2), 2 => (1, 1, 2), 1 => (1, 1, 1))
# core_config = DataStructures.OrderedDict(16 => (2, 2, 4), 12 => (2, 2, 3), 8 => (2, 2, 2), 4 => (1, 2, 2), 2 => (1, 1, 2), 1 => (1, 1, 1))
node_core_partitions = DataStructures.OrderedDict((node, core) => (node_partition .* core_partition) for (core, core_partition) in core_config for (node, node_partition) in node_config)
