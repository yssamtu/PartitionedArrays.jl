include("../test_fdm.jl")
nparts = (2,2,2)
distributed_run(test_fdm,mpi,nparts)

