module MatrixAssemblyBenchmarks

using PartitionedArrays
using PartitionedArrays: local_permutation, split_matrix_blocks, split_matrix
using SparseArrays
using PetscCall
import JSON
using Base64

include("../../src/matrix_assembly.jl")
include("benchmarks.jl")

end # module