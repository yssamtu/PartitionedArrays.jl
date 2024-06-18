module MatrixAssemblyBenchmarks

using PartitionedArrays
using PartitionedArrays: local_permutation, split_matrix_blocks, split_matrix
using SparseArrays
using PetscCall
using MPI
import DataStructures
import JSON

include("../../src/matrix_assembly.jl")
include("common.jl")
include("helper.jl")
include("benchmarks.jl")

end # module