function assemble_matrix_no_compressed_snd_and_no_cache!(I, J, V, rows, cols, f)
    function quick_sort_partition!(part, key, values::Vararg{Any,N}) where {N}
        global_to_own_part = global_to_own(part)
        left_ptr = firstindex(key)
        right_ptr = lastindex(key)
        while true
            while left_ptr <= right_ptr && global_to_own_part[key[left_ptr]] != 0
                left_ptr += 1
            end
            while left_ptr <= right_ptr && global_to_own_part[key[right_ptr]] == 0
                right_ptr -= 1
            end
            if left_ptr > right_ptr
                break
            end
            key[left_ptr], key[right_ptr] = key[right_ptr], key[left_ptr]
            for i in 1:N
                values[i][left_ptr], values[i][right_ptr] = values[i][right_ptr], values[i][left_ptr]
            end
            left_ptr += 1
            right_ptr -= 1
        end
        right_ptr
    end
    function partition_and_setup_cache_snd!(I, J, V, I_owner, parts_snd, rows_sa)
        n_hold_data = quick_sort_partition!(rows_sa, I, J, V, I_owner)
        snd_index = (n_hold_data+1):lastindex(I)
        I_raw_snd_data = view(I, snd_index)
        J_raw_snd_data = view(J, snd_index)
        V_raw_snd_data = view(V, snd_index)
        I_raw_snd_owner = view(I_owner, snd_index)
        n_snd_data = length(I) - n_hold_data
        I_snd_data = similar(I, n_snd_data)
        J_snd_data = similar(I, n_snd_data)
        V_snd_data = similar(V, n_snd_data)
        owner_to_p = Dict(owner => i for (i, owner) in enumerate(parts_snd))
        ptrs = zeros(Int32, length(parts_snd) + 1)
        for (i, owner) in enumerate(I_raw_snd_owner)
            p = owner_to_p[owner]
            I_raw_snd_owner[i] = p
            ptrs[p+1] += 1
        end
        length_to_ptrs!(ptrs)
        for (i, j, v, p) in zip(I_raw_snd_data, J_raw_snd_data, V_raw_snd_data, I_raw_snd_owner)
            index = ptrs[p]
            I_snd_data[index] = i
            J_snd_data[index] = j
            V_snd_data[index] = v
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data
    end
    function store_recv_data!(I, J, V, n_hold_data, I_rcv, J_rcv, V_rcv)
        n_data = n_hold_data + length(I_rcv.data)
        resize!(I, n_data)
        sizehint!(I, n_data)
        resize!(J, n_data)
        sizehint!(J, n_data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        I[rcv_index] = I_rcv.data
        J[rcv_index] = J_rcv.data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(I, J, V, rows_fa, cols_fa)
        n_own_data = quick_sort_partition!(cols_fa, J, I, V)
        is_own = firstindex(I):n_own_data
        is_ghost = (n_own_data+1):lastindex(I)
        I_own_own = view(I, is_own)
        J_own_own = view(J, is_own)
        V_own_own = view(V, is_own)
        I_own_ghost = view(I, is_ghost)
        J_own_ghost = view(J, is_ghost)
        V_own_ghost = view(V, is_ghost)
        map_global_to_own!(I_own_own, rows_fa)
        map_global_to_own!(J_own_own, cols_fa)
        map_global_to_own!(I_own_ghost, rows_fa)
        map_global_to_ghost!(J_own_ghost, cols_fa)
        n_own_rows = own_length(rows_fa)
        n_own_cols = own_length(cols_fa)
        n_ghost_rows = ghost_length(rows_fa)
        n_ghost_cols = ghost_length(cols_fa)
        Ti = eltype(I)
        Tv = eltype(V)
        combine = +
        own_own = f(I_own_own, J_own_own, V_own_own, n_own_rows, n_own_cols, combine)
        own_ghost = f(I_own_ghost, J_own_ghost, V_own_ghost, n_own_rows, n_ghost_cols, combine)
        ghost_own = f(Ti[], Ti[], Tv[], n_ghost_rows, n_own_cols, combine)
        ghost_ghost = f(Ti[], Ti[], Tv[], n_ghost_rows, n_ghost_cols, combine)
        blocks = split_matrix_blocks(own_own, own_ghost, ghost_own, ghost_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm)
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size = map(partition_and_setup_cache_snd!, I, J, V, I_owner, parts_snd, rows_sa) |> tuple_of_arrays
    graph = ExchangeGraph(parts_snd, parts_rcv)
    t_I = exchange(I_snd_buf, graph)
    t_J = exchange(J_snd_buf, graph)
    t_V = exchange(V_snd_buf, graph)
    @async begin
        I_rcv_buf = fetch(t_I)
        J_rcv_buf = fetch(t_J)
        V_rcv_buf = fetch(t_V)
        map(store_recv_data!, I, J, V, hold_data_size, I_rcv_buf, J_rcv_buf, V_rcv_buf)
        rows_fa = rows
        J_owner = find_owner(cols, J)
        cols_fa = map(union_ghost, cols, J, J_owner)
        vals_fa = map(split_and_compress!, I, J, V, rows_fa, cols_fa)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled)
    end
end

function assemble_matrix_with_compressed_snd_and_no_cache!(I, J, V, rows, cols, f)
    function quick_sort_partition!(part, key, values::Vararg{Any,N}) where {N}
        global_to_own_part = global_to_own(part)
        left_ptr = firstindex(key)
        right_ptr = lastindex(key)
        while true
            while left_ptr <= right_ptr && global_to_own_part[key[left_ptr]] != 0
                left_ptr += 1
            end
            while left_ptr <= right_ptr && global_to_own_part[key[right_ptr]] == 0
                right_ptr -= 1
            end
            if left_ptr > right_ptr
                break
            end
            key[left_ptr], key[right_ptr] = key[right_ptr], key[left_ptr]
            for i in 1:N
                values[i][left_ptr], values[i][right_ptr] = values[i][right_ptr], values[i][left_ptr]
            end
            left_ptr += 1
            right_ptr -= 1
        end
        right_ptr
    end
    function partition_and_setup_cache_snd!(I, J, V, parts_snd, rows_sa, cols)
        n_hold_data = quick_sort_partition!(rows_sa, I, J, V)
        snd_index = (n_hold_data+1):lastindex(I)
        I_raw_snd_data = view(I, snd_index)
        J_raw_snd_data = view(J, snd_index)
        V_raw_snd_data = view(V, snd_index)
        map_global_to_ghost!(I_raw_snd_data, rows_sa)
        n_ghost_rows = ghost_length(rows_sa)
        n_global_cols = global_length(cols)
        compressed_snd_data = sparse(J_raw_snd_data, I_raw_snd_data, V_raw_snd_data, n_global_cols, n_ghost_rows)
        n_snd_data = nnz(compressed_snd_data)
        I_snd_data = similar(I, n_snd_data)
        J_snd_data = similar(J, n_snd_data)
        V_snd_data = similar(V, n_snd_data)
        owner_to_p = Dict(owner => i for (i, owner) in enumerate(parts_snd))
        ghost_to_p = [owner_to_p[owner] for owner in ghost_to_owner(rows_sa)]
        ptrs = zeros(Int32, length(parts_snd) + 1)
        for (_, i, _) in nziterator(compressed_snd_data)
            ptrs[ghost_to_p[i]+1] += 1
        end
        length_to_ptrs!(ptrs)
        ghost_to_global_row = ghost_to_global(rows_sa)
        for (j, i, v) in nziterator(compressed_snd_data)
            p = ghost_to_p[i]
            index = ptrs[p]
            I_snd_data[index] = ghost_to_global_row[i]
            J_snd_data[index] = j
            V_snd_data[index] = v
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data
    end
    function store_recv_data!(I, J, V, n_hold_data, I_rcv, J_rcv, V_rcv)
        n_data = n_hold_data + length(I_rcv.data)
        resize!(I, n_data)
        sizehint!(I, n_data)
        resize!(J, n_data)
        sizehint!(J, n_data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        I[rcv_index] = I_rcv.data
        J[rcv_index] = J_rcv.data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(I, J, V, rows_fa, cols_fa)
        n_own_data = quick_sort_partition!(cols_fa, J, I, V)
        is_own = firstindex(I):n_own_data
        is_ghost = (n_own_data+1):lastindex(I)
        I_own_own = view(I, is_own)
        J_own_own = view(J, is_own)
        V_own_own = view(V, is_own)
        I_own_ghost = view(I, is_ghost)
        J_own_ghost = view(J, is_ghost)
        V_own_ghost = view(V, is_ghost)
        map_global_to_own!(I_own_own, rows_fa)
        map_global_to_own!(J_own_own, cols_fa)
        map_global_to_own!(I_own_ghost, rows_fa)
        map_global_to_ghost!(J_own_ghost, cols_fa)
        n_own_rows = own_length(rows_fa)
        n_own_cols = own_length(cols_fa)
        n_ghost_rows = ghost_length(rows_fa)
        n_ghost_cols = ghost_length(cols_fa)
        Ti = eltype(I)
        Tv = eltype(V)
        combine = +
        own_own = f(I_own_own, J_own_own, V_own_own, n_own_rows, n_own_cols, combine)
        own_ghost = f(I_own_ghost, J_own_ghost, V_own_ghost, n_own_rows, n_ghost_cols, combine)
        ghost_own = f(Ti[], Ti[], Tv[], n_ghost_rows, n_own_cols, combine)
        ghost_ghost = f(Ti[], Ti[], Tv[], n_ghost_rows, n_ghost_cols, combine)
        blocks = split_matrix_blocks(own_own, own_ghost, ghost_own, ghost_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm)
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size = map(partition_and_setup_cache_snd!, I, J, V, parts_snd, rows_sa, cols) |> tuple_of_arrays
    graph = ExchangeGraph(parts_snd, parts_rcv)
    t_I = exchange(I_snd_buf, graph)
    t_J = exchange(J_snd_buf, graph)
    t_V = exchange(V_snd_buf, graph)
    @async begin
        I_rcv_buf = fetch(t_I)
        J_rcv_buf = fetch(t_J)
        V_rcv_buf = fetch(t_V)
        map(store_recv_data!, I, J, V, hold_data_size, I_rcv_buf, J_rcv_buf, V_rcv_buf)
        rows_fa = rows
        J_owner = find_owner(cols, J)
        cols_fa = map(union_ghost, cols, J, J_owner)
        vals_fa = map(split_and_compress!, I, J, V, rows_fa, cols_fa)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled)
    end
end