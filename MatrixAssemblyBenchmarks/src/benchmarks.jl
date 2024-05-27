function coo_scalar_fem(cells_per_dir, parts_per_dir, parts, ::Type{T}, ::Type{Ti}) where {Ti,T}
    # TODO only symbolic info for the moment
    D = length(cells_per_dir)
    nodes_per_dir = cells_per_dir .+ 1
    ghost_per_dir = ntuple(d -> true, Val(D))
    node_partition = uniform_partition(parts, parts_per_dir, nodes_per_dir, ghost_per_dir)
    cell_partition = uniform_partition(parts, parts_per_dir, cells_per_dir)
    isboundary(cartesian_node) = any(map((d, i) -> i ∉ 2:(nodes_per_dir[d]-1), 1:D, Tuple(cartesian_node)))
    function fill_mask!(local_node_to_mask, nodes)
        node_to_cartesian_node = CartesianIndices(nodes_per_dir)
        n_local_nodes = local_length(nodes)
        local_node_to_node = local_to_global(nodes)
        for local_node in 1:n_local_nodes
            node = local_node_to_node[local_node]
            cartesian_node = node_to_cartesian_node[node]
            mask = !isboundary(cartesian_node)
            local_node_to_mask[local_node] = mask
        end
    end
    node_to_mask = pfill(false, node_partition)
    map(fill_mask!, partition(node_to_mask), node_partition)
    dof_to_local_node, node_to_local_dof = find_local_indices(node_to_mask)
    dof_partition = partition(axes(dof_to_local_node, 1))
    function setup(cells, nodes, dofs, local_node_to_local_dof)
        cell_to_cartesian_cell = CartesianIndices(cells_per_dir)
        cartesian_node_to_node = LinearIndices(nodes_per_dir)
        own_cell_to_cell = own_to_global(cells)
        node_to_local_node = global_to_local(nodes)
        local_dof_to_dof = local_to_global(dofs)
        n_own_cells = length(own_cell_to_cell)
        n_nz = n_own_cells * (2^(2 * D))
        myI = zeros(Ti, n_nz)
        myJ = zeros(Ti, n_nz)
        myV = ones(T, n_nz)
        p = 0
        for own_cell in 1:n_own_cells
            cell = own_cell_to_cell[own_cell]
            cartesian_cell = cell_to_cartesian_cell[cell]
            for di in 1:D
                offset_i = ntuple(d -> (d == di ? 1 : 0), Val(D))
                cartesian_node_i = CartesianIndex(Tuple(cartesian_cell) .+ offset_i)
                if isboundary(cartesian_node_i)
                    continue
                end
                node_i = cartesian_node_to_node[cartesian_node_i]
                local_node_i = node_to_local_node[node_i]
                local_dof_i = local_node_to_local_dof[local_node_i]
                dof_i = local_dof_to_dof[local_dof_i]
                for dj in 1:D
                    offset_j = ntuple(d -> (d == dj ? 1 : 0), Val(D))
                    cartesian_node_j = CartesianIndex(Tuple(cartesian_cell) .+ offset_j)
                    if isboundary(cartesian_node_j)
                        continue
                    end
                    node_j = cartesian_node_to_node[cartesian_node_j]
                    local_node_j = node_to_local_node[node_j]
                    local_dof_j = local_node_to_local_dof[local_node_j]
                    dof_j = local_dof_to_dof[local_dof_j]
                    p += 1
                    myI[p] = dof_i
                    myJ[p] = dof_j
                end
            end
        end
        (myI[1:p], myJ[1:p], myV[1:p])
    end
    I, J, V = map(setup, cell_partition, node_partition, dof_partition, partition(node_to_local_dof)) |> tuple_of_arrays
    row_partition = map(remove_ghost, dof_partition)
    col_partition = row_partition
    I, J, V, row_partition, col_partition
end

function benchmark_psparse(distribute, params)
    nruns, cells_per_dir, parts_per_dir, method = params
    np = prod(parts_per_dir)
    parts = distribute(LinearIndices((np,)))
    t_buildmat = zeros(nruns)
    t_rebuildmat = zeros(nruns)
    petsc_comm = PetscCall.setup_petsc_comm(parts)
    function petsc_setvalues(I, J, V, rows, cols)
        m = own_length(rows)
        n = own_length(cols)
        M = global_length(rows)
        N = global_length(cols)
        I .= I .- 1
        J .= J .- 1
        for irun in 1:nruns
            t_buildmat[irun] = @elapsed begin
                A = Ref{PetscCall.Mat}()
                PetscCall.@check_error_code PetscCall.MatCreate(petsc_comm, A)
                PetscCall.@check_error_code PetscCall.MatSetType(A[], PetscCall.MATMPIAIJ)
                PetscCall.@check_error_code PetscCall.MatSetSizes(A[], m, n, M, N)
                PetscCall.@check_error_code PetscCall.MatMPIAIJSetPreallocation(A[], 32, C_NULL, 32, C_NULL)
                for p in 1:length(I)
                    PetscCall.@check_error_code PetscCall.MatSetValues(A[], 1, view(I, p:p), 1, view(J, p:p), view(V, p:p), PetscCall.ADD_VALUES)
                end
                PetscCall.@check_error_code PetscCall.MatAssemblyBegin(A[], PetscCall.MAT_FINAL_ASSEMBLY)
                PetscCall.@check_error_code PetscCall.MatAssemblyEnd(A[], PetscCall.MAT_FINAL_ASSEMBLY)
            end
            PetscCall.@check_error_code PetscCall.MatDestroy(A)
        end
        A = Ref{PetscCall.Mat}()
        PetscCall.@check_error_code PetscCall.MatCreate(petsc_comm, A)
        PetscCall.@check_error_code PetscCall.MatSetType(A[], PetscCall.MATMPIAIJ)
        PetscCall.@check_error_code PetscCall.MatSetSizes(A[], m, n, M, N)
        PetscCall.@check_error_code PetscCall.MatMPIAIJSetPreallocation(A[], 32, C_NULL, 32, C_NULL)
        for p in 1:length(I)
            PetscCall.@check_error_code PetscCall.MatSetValues(A[], 1, view(I, p:p), 1, view(J, p:p), view(V, p:p), PetscCall.ADD_VALUES)
        end
        PetscCall.@check_error_code PetscCall.MatAssemblyBegin(A[], PetscCall.MAT_FINAL_ASSEMBLY)
        PetscCall.@check_error_code PetscCall.MatAssemblyEnd(A[], PetscCall.MAT_FINAL_ASSEMBLY)
        for irun in 1:nruns
            t_rebuildmat[irun] = @elapsed begin
                for p in 1:length(I)
                    PetscCall.@check_error_code PetscCall.MatSetValues(A[], 1, view(I, p:p), 1, view(J, p:p), view(V, p:p), PetscCall.ADD_VALUES)
                end
                PetscCall.@check_error_code PetscCall.MatAssemblyBegin(A[], PetscCall.MAT_FINAL_ASSEMBLY)
                PetscCall.@check_error_code PetscCall.MatAssemblyEnd(A[], PetscCall.MAT_FINAL_ASSEMBLY)
            end
        end
        PetscCall.@check_error_code PetscCall.MatDestroy(A)
    end
    function petsc_coo(I, J, V, rows, cols)
        m = own_length(rows)
        n = own_length(cols)
        M = global_length(rows)
        N = global_length(cols)
        I .= I .- 1
        J .= J .- 1
        ncoo = length(I)
        for irun in 1:nruns
            t_buildmat[irun] = @elapsed begin
                A = Ref{PetscCall.Mat}()
                PetscCall.@check_error_code PetscCall.MatCreate(petsc_comm, A)
                PetscCall.@check_error_code PetscCall.MatSetType(A[], PetscCall.MATMPIAIJ)
                PetscCall.@check_error_code PetscCall.MatSetSizes(A[], m, n, M, N)
                PetscCall.@check_error_code PetscCall.MatSetPreallocationCOO(A[], ncoo, I, J)
                PetscCall.@check_error_code PetscCall.MatSetValuesCOO(A[], V, PetscCall.ADD_VALUES)
                PetscCall.@check_error_code PetscCall.MatAssemblyBegin(A[], PetscCall.MAT_FINAL_ASSEMBLY)
                PetscCall.@check_error_code PetscCall.MatAssemblyEnd(A[], PetscCall.MAT_FINAL_ASSEMBLY)
            end
            PetscCall.@check_error_code PetscCall.MatDestroy(A)
        end
        A = Ref{PetscCall.Mat}()
        PetscCall.@check_error_code PetscCall.MatCreate(petsc_comm, A)
        PetscCall.@check_error_code PetscCall.MatSetType(A[], PetscCall.MATMPIAIJ)
        PetscCall.@check_error_code PetscCall.MatSetSizes(A[], m, n, M, N)
        PetscCall.@check_error_code PetscCall.MatSetPreallocationCOO(A[], ncoo, I, J)
        PetscCall.@check_error_code PetscCall.MatSetValuesCOO(A[], V, PetscCall.ADD_VALUES)
        PetscCall.@check_error_code PetscCall.MatAssemblyBegin(A[], PetscCall.MAT_FINAL_ASSEMBLY)
        PetscCall.@check_error_code PetscCall.MatAssemblyEnd(A[], PetscCall.MAT_FINAL_ASSEMBLY)
        for irun in 1:nruns
            t_rebuildmat[irun] = @elapsed begin
                PetscCall.@check_error_code PetscCall.MatSetValuesCOO(A[], V, PetscCall.ADD_VALUES)
                PetscCall.@check_error_code PetscCall.MatAssemblyBegin(A[], PetscCall.MAT_FINAL_ASSEMBLY)
                PetscCall.@check_error_code PetscCall.MatAssemblyEnd(A[], PetscCall.MAT_FINAL_ASSEMBLY)
            end
        end
        PetscCall.@check_error_code PetscCall.MatDestroy(A)
    end
    Ti = PetscCall.PetscInt
    T = PetscCall.PetscScalar
    psparse_args = coo_scalar_fem(cells_per_dir, parts_per_dir, parts, T, Ti)
    V = psparse_args[3]
    if method == "psparse"
        for irun in 1:nruns
            t_buildmat[irun] = @elapsed psparse(psparse_args...; reuse=true) |> fetch
        end
        A, cacheA = psparse(psparse_args...; reuse=true) |> fetch
        for irun in 1:nruns
            t_rebuildmat[irun] = @elapsed psparse!(A, V, cacheA) |> wait
        end
    elseif method == "petsc_setvalues"
        PetscCall.init(finalize_atexit=false)
        map(petsc_setvalues, psparse_args...)
        PetscCall.finalize()
    elseif method == "petsc_coo"
        PetscCall.init(finalize_atexit=false)
        map(petsc_coo, psparse_args...)
        PetscCall.finalize()
    elseif method == "assemble_matrix_no_compressed_snd_and_with_int_vector_cache"
        for irun in 1:nruns
            copy_psparse_args = deepcopy(psparse_args)
            t_buildmat[irun] = @elapsed assemble_matrix_no_compressed_snd_and_with_int_vector_cache!(sparse, copy_psparse_args...) |> fetch
        end
        copy_psparse_args = deepcopy(psparse_args)
        A, cacheA = assemble_matrix_no_compressed_snd_and_with_int_vector_cache!(sparse, copy_psparse_args...) |> fetch
        for irun in 1:nruns
            copy_V = copy(psparse_args[3])
            t_rebuildmat[irun] = @elapsed assemble_matrix_no_compressed_snd_and_with_int_vector_cache!(A, copy_V, cacheA) |> wait
        end
    elseif method == "assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache"
        for irun in 1:nruns
            copy_psparse_args = deepcopy(psparse_args)
            t_buildmat[irun] = @elapsed assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache!(sparse, copy_psparse_args...) |> fetch
        end
        copy_psparse_args = deepcopy(psparse_args)
        A, cacheA = assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache!(sparse, copy_psparse_args...) |> fetch
        for irun in 1:nruns
            copy_V = copy(psparse_args[3])
            t_rebuildmat[irun] = @elapsed assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache!(A, copy_V, cacheA) |> wait
        end
    elseif method == "assemble_matrix_no_compressed_snd_and_with_auto_cache"
        for irun in 1:nruns
            copy_psparse_args = deepcopy(psparse_args)
            t_buildmat[irun] = @elapsed assemble_matrix_no_compressed_snd_and_with_auto_cache!(sparse, copy_psparse_args...) |> fetch
        end
        copy_psparse_args = deepcopy(psparse_args)
        A, cacheA = assemble_matrix_no_compressed_snd_and_with_auto_cache!(sparse, copy_psparse_args...) |> fetch
        for irun in 1:nruns
            copy_V = copy(psparse_args[3])
            t_rebuildmat[irun] = @elapsed assemble_matrix_no_compressed_snd_and_with_auto_cache!(A, copy_V, cacheA) |> wait
        end
    elseif method == "assemble_matrix_with_compressed_snd_and_with_int_vector_cache"
        for irun in 1:nruns
            copy_psparse_args = deepcopy(psparse_args)
            t_buildmat[irun] = @elapsed assemble_matrix_with_compressed_snd_and_with_int_vector_cache!(sparse, copy_psparse_args...) |> fetch
        end
        copy_psparse_args = deepcopy(psparse_args)
        A, cacheA = assemble_matrix_with_compressed_snd_and_with_int_vector_cache!(sparse, copy_psparse_args...) |> fetch
        for irun in 1:nruns
            copy_V = copy(psparse_args[3])
            t_rebuildmat[irun] = @elapsed assemble_matrix_with_compressed_snd_and_with_int_vector_cache!(A, copy_V, cacheA) |> wait
        end
    elseif method == "assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache"
        for irun in 1:nruns
            copy_psparse_args = deepcopy(psparse_args)
            t_buildmat[irun] = @elapsed assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache!(sparse, copy_psparse_args...) |> fetch
        end
        copy_psparse_args = deepcopy(psparse_args)
        A, cacheA = assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache!(sparse, copy_psparse_args...) |> fetch
        for irun in 1:nruns
            copy_V = copy(psparse_args[3])
            t_rebuildmat[irun] = @elapsed assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache!(A, copy_V, cacheA) |> wait
        end
    elseif method == "assemble_matrix_with_compressed_snd_and_with_auto_cache"
        for irun in 1:nruns
            copy_psparse_args = deepcopy(psparse_args)
            t_buildmat[irun] = @elapsed assemble_matrix_with_compressed_snd_and_with_auto_cache!(sparse, copy_psparse_args...) |> fetch
        end
        copy_psparse_args = deepcopy(psparse_args)
        A, cacheA = assemble_matrix_with_compressed_snd_and_with_auto_cache!(sparse, copy_psparse_args...) |> fetch
        for irun in 1:nruns
            copy_V = copy(psparse_args[3])
            t_rebuildmat[irun] = @elapsed assemble_matrix_with_compressed_snd_and_with_auto_cache!(A, copy_V, cacheA) |> wait
        end
    end
    ts_in_main = gather(map(p -> (; t_buildmat, t_rebuildmat), parts))
    results_in_main = map_main(ts_in_main) do ts
        buildmat = map(i -> i.t_buildmat, ts)
        rebuildmat = map(i -> i.t_rebuildmat, ts)
        results = (; params..., buildmat, rebuildmat)
        results
    end
end

function get_job_name(params)
    params_value = values(params)
    params_str = sprint(show, params_value)
    params_str = replace(params_str, " " => "")
    base64_str = base64encode(params_str)
    replace(base64_str, "+" => "-", "/" => "_")
end

function get_job_name(file_name::String)
    file_base_name = basename(file_name)
    splitext(file_base_name)[1]
end

function get_file_name(params)
    folder = mkpath("results")
    extension = ".json"
    job_name = get_job_name(params)
    joinpath(folder, job_name * extension)
end

function get_params(file_name::String)
    function parse_tuple(type, str)
        vec_str = split(str, ",")
        vec_type = parse.(type, vec_str)
        NTuple{length(vec_type),eltype(vec_type)}(vec_type)
    end
    job_name = get_job_name(file_name)
    base64_str = replace(job_name, "-" => "+", "_" => "/")
    params_str = String(base64decode(base64_str))
    params_list = split(replace(params_str, ")" => "", "\"" => "(")[2:end-1], ",(")
    nruns = parse(Int, params_list[1])
    cells_per_dir = parse_tuple(Int, params_list[2])
    parts_per_dir = parse_tuple(Int, params_list[3])
    method = String(params_list[4])
    (; nruns, cells_per_dir, parts_per_dir, method)
end

function experiment(params)
    results_in_main = with_mpi(distribute -> benchmark_psparse(distribute, params))
    map_main(results_in_main) do results
        open(get_file_name(params), "w") do f
            JSON.print(f, results, 2)
        end
    end
end