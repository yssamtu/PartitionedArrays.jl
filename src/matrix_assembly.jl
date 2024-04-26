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

function precompute_nzindex!(K::AbstractVector{Int32}, A, I, J)
    for (p, (i, j)) in enumerate(zip(I, J))
        if i < 1 || j < 1
            continue
        end
        K[p] = nzindex(A, i, j)
    end
end

function assemble_matrix_no_compressed_snd_and_with_int_vector_cache!(I, J, V, rows, cols, f)
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
            left_ptr += 1
            right_ptr -= 1
        end
        Tkey = eltype(key)
        change = collect(Tkey, 1:right_ptr)
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
            change[left_ptr] = right_ptr
            left_ptr += 1
            right_ptr -= 1
        end
        right_ptr, change
    end
    function partition_and_setup_cache_snd!(I, J, V, I_owner, parts_snd, rows_sa)
        n_hold_data, change_index = quick_sort_partition!(rows_sa, I, J, V, I_owner)
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
        for (n, (i, j, v, p)) in enumerate(zip(I_raw_snd_data, J_raw_snd_data, V_raw_snd_data, I_raw_snd_owner))
            index = ptrs[p]
            I_snd_data[index] = i
            J_snd_data[index] = j
            V_snd_data[index] = v
            I_owner[n] = index
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        resize!(I_owner, n_snd_data)
        sizehint!(I_owner, n_snd_data)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data, change_index, I_owner
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
    function split_and_compress!(I, J, V, perm, rows_fa, cols_fa)
        n_own_data, change_index = quick_sort_partition!(cols_fa, J, I, V)
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
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        precompute_nzindex!(perm_own, own_own, I_own_own, J_own_own)
        precompute_nzindex!(perm_ghost, own_ghost, I_own_ghost, J_own_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm), change_index, perm
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size, change_snd, perm_snd = map(partition_and_setup_cache_snd!, I, J, V, I_owner, parts_snd, rows_sa) |> tuple_of_arrays
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
        vals_fa, change_sparse, perm_sparse = map(split_and_compress!, I, J, V, J_owner, rows_fa, cols_fa) |> tuple_of_arrays
        cache = (; graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled), cache
    end
end

function assemble_matrix_no_compressed_snd_and_with_int_vector_cache!(A, V, cache)
    function perm_partition!(V, perm::Vector{T}) where {T}
        for (i, p) in enumerate(perm)
            V[i], V[p] = V[p], V[i]
        end
    end
    function partition_and_setup_cache_snd!(V_snd, V, n_hold_data, change_index, perm)
        perm_partition!(V, change_index)
        snd_index = (n_hold_data+1):lastindex(V)
        V_raw_snd_data = view(V, snd_index)
        V_snd_data = V_snd.data
        for (p, v) in zip(perm, V_raw_snd_data)
            V_snd_data[p] = v
        end
    end
    function store_recv_data!(V, n_hold_data, V_rcv)
        n_data = n_hold_data + length(V_rcv.data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(A, V, change_index, perm)
        perm_partition!(V, change_index)
        n_own_data = length(change_index)
        is_own = firstindex(V):n_own_data
        is_ghost = (n_own_data+1):lastindex(V)
        V_own_own = view(V, is_own)
        V_own_ghost = view(V, is_ghost)
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        sparse_matrix!(A.blocks.own_own, V_own_own, perm_own; reset=false)
        sparse_matrix!(A.blocks.own_ghost, V_own_ghost, perm_ghost; reset=false)
        return
    end
    graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse = cache
    map(partition_and_setup_cache_snd!, V_snd_buf, V, hold_data_size, change_snd, perm_snd)
    t_V = exchange!(V_rcv_buf, V_snd_buf, graph)
    @async begin
        fetch(t_V)
        map(store_recv_data!, V, hold_data_size, V_rcv_buf)
        map(split_and_compress!, partition(A), V, change_sparse, perm_sparse)
        A
    end
end

function assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache!(I, J, V, rows, cols, f)
    function quick_sort_partition!(part, key, values::Vararg{Any,N}) where {N}
        global_to_own_part = global_to_own(part)
        n_change = 0
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
            n_change += 1
            left_ptr += 1
            right_ptr -= 1
        end
        Tkey = eltype(key)
        change = Vector{Tuple{Tkey,Tkey}}(undef, n_change)
        ptr = firstindex(change)
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
            change[ptr] = (left_ptr, right_ptr)
            ptr += 1
            left_ptr += 1
            right_ptr -= 1
        end
        right_ptr, change
    end
    function partition_and_setup_cache_snd!(I, J, V, I_owner, parts_snd, rows_sa)
        n_hold_data, change_index = quick_sort_partition!(rows_sa, I, J, V, I_owner)
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
        for (n, (i, j, v, p)) in enumerate(zip(I_raw_snd_data, J_raw_snd_data, V_raw_snd_data, I_raw_snd_owner))
            index = ptrs[p]
            I_snd_data[index] = i
            J_snd_data[index] = j
            V_snd_data[index] = v
            I_owner[n] = index
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        resize!(I_owner, n_snd_data)
        sizehint!(I_owner, n_snd_data)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data, change_index, I_owner
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
    function split_and_compress!(I, J, V, perm, rows_fa, cols_fa)
        n_own_data, change_index = quick_sort_partition!(cols_fa, J, I, V)
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
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        precompute_nzindex!(perm_own, own_own, I_own_own, J_own_own)
        precompute_nzindex!(perm_ghost, own_ghost, I_own_ghost, J_own_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm), change_index, perm
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size, change_snd, perm_snd = map(partition_and_setup_cache_snd!, I, J, V, I_owner, parts_snd, rows_sa) |> tuple_of_arrays
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
        vals_fa, change_sparse, perm_sparse = map(split_and_compress!, I, J, V, J_owner, rows_fa, cols_fa) |> tuple_of_arrays
        cache = (; graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled), cache
    end
end

function assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache!(A, V, cache)
    function perm_partition!(V, perm::Vector{Tuple{T,T}}) where {T}
        for (i, j) in perm
            V[i], V[j] = V[j], V[i]
        end
    end
    function partition_and_setup_cache_snd!(V_snd, V, n_hold_data, change_index, perm)
        perm_partition!(V, change_index)
        snd_index = (n_hold_data+1):lastindex(V)
        V_raw_snd_data = view(V, snd_index)
        V_snd_data = V_snd.data
        for (p, v) in zip(perm, V_raw_snd_data)
            V_snd_data[p] = v
        end
    end
    function store_recv_data!(V, n_hold_data, V_rcv)
        n_data = n_hold_data + length(V_rcv.data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(A, V, change_index, perm)
        perm_partition!(V, change_index)
        n_own_data = isempty(change_index) ? 1 : last(change_index)[1]
        is_own = firstindex(V):n_own_data
        is_ghost = (n_own_data+1):lastindex(V)
        V_own_own = view(V, is_own)
        V_own_ghost = view(V, is_ghost)
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        sparse_matrix!(A.blocks.own_own, V_own_own, perm_own; reset=false)
        sparse_matrix!(A.blocks.own_ghost, V_own_ghost, perm_ghost; reset=false)
        return
    end
    graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse = cache
    map(partition_and_setup_cache_snd!, V_snd_buf, V, hold_data_size, change_snd, perm_snd)
    t_V = exchange!(V_rcv_buf, V_snd_buf, graph)
    @async begin
        fetch(t_V)
        map(store_recv_data!, V, hold_data_size, V_rcv_buf)
        map(split_and_compress!, partition(A), V, change_sparse, perm_sparse)
        A
    end
end

function assemble_matrix_no_compressed_snd_and_with_auto_cache!(I, J, V, rows, cols, f)
    function quick_sort_partition!(part, key, values::Vararg{Any,N}) where {N}
        global_to_own_part = global_to_own(part)
        n_change = 0
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
            n_change += 1
            left_ptr += 1
            right_ptr -= 1
        end
        Tkey = eltype(key)
        if right_ptr - n_change - n_change <= 0
            change = collect(Tkey, 1:right_ptr)
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
                change[left_ptr] = right_ptr
                left_ptr += 1
                right_ptr -= 1
            end
        else
            change = Vector{Tuple{Tkey,Tkey}}(undef, n_change)
            ptr = firstindex(change)
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
                change[ptr] = (left_ptr, right_ptr)
                ptr += 1
                left_ptr += 1
                right_ptr -= 1
            end
        end
        right_ptr, change
    end
    function partition_and_setup_cache_snd!(I, J, V, I_owner, parts_snd, rows_sa)
        n_hold_data, change_index = quick_sort_partition!(rows_sa, I, J, V, I_owner)
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
        for (n, (i, j, v, p)) in enumerate(zip(I_raw_snd_data, J_raw_snd_data, V_raw_snd_data, I_raw_snd_owner))
            index = ptrs[p]
            I_snd_data[index] = i
            J_snd_data[index] = j
            V_snd_data[index] = v
            I_owner[n] = index
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        resize!(I_owner, n_snd_data)
        sizehint!(I_owner, n_snd_data)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data, change_index, I_owner
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
    function split_and_compress!(I, J, V, perm, rows_fa, cols_fa)
        n_own_data, change_index = quick_sort_partition!(cols_fa, J, I, V)
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
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        precompute_nzindex!(perm_own, own_own, I_own_own, J_own_own)
        precompute_nzindex!(perm_ghost, own_ghost, I_own_ghost, J_own_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm), change_index, perm
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size, change_snd, perm_snd = map(partition_and_setup_cache_snd!, I, J, V, I_owner, parts_snd, rows_sa) |> tuple_of_arrays
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
        vals_fa, change_sparse, perm_sparse = map(split_and_compress!, I, J, V, J_owner, rows_fa, cols_fa) |> tuple_of_arrays
        cache = (; graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled), cache
    end
end

function assemble_matrix_no_compressed_snd_and_with_auto_cache!(A, V, cache)
    function partition_and_setup_cache_snd!(V_snd, V, n_hold_data, change_index, perm)
        perm_partition!(v, change::Vector{T}) where {T} = begin
            for (i, p) in enumerate(change)
                v[i], v[p] = v[p], v[i]
            end
        end
        perm_partition!(v, change::Vector{Tuple{T,T}}) where {T} = begin
            for (i, j) in change
                v[i], v[j] = v[j], v[i]
            end
        end
        perm_partition!(V, change_index)
        snd_index = (n_hold_data+1):lastindex(V)
        V_raw_snd_data = view(V, snd_index)
        V_snd_data = V_snd.data
        for (p, v) in zip(perm, V_raw_snd_data)
            V_snd_data[p] = v
        end
    end
    function store_recv_data!(V, n_hold_data, V_rcv)
        n_data = n_hold_data + length(V_rcv.data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(A, V, change_index, perm)
        perm_partition!(V, perm::Vector{T}) where {T} = begin
            for (i, p) in enumerate(perm)
                V[i], V[p] = V[p], V[i]
            end
        end
        perm_partition!(V, perm::Vector{Tuple{T,T}}) where {T} = begin
            for (i, j) in perm
                V[i], V[j] = V[j], V[i]
            end
        end
        get_n_own_data(change_index::Vector{T}) where {T} = length(change_index)
        get_n_own_data(change_index::Vector{Tuple{T,T}}) where {T} = isempty(change_index) ? 1 : last(change_index)[1]
        perm_partition!(V, change_index)
        n_own_data = get_n_own_data(change_index)
        is_own = firstindex(V):n_own_data
        is_ghost = (n_own_data+1):lastindex(V)
        V_own_own = view(V, is_own)
        V_own_ghost = view(V, is_ghost)
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        sparse_matrix!(A.blocks.own_own, V_own_own, perm_own; reset=false)
        sparse_matrix!(A.blocks.own_ghost, V_own_ghost, perm_ghost; reset=false)
        return
    end
    graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse = cache
    map(partition_and_setup_cache_snd!, V_snd_buf, V, hold_data_size, change_snd, perm_snd)
    t_V = exchange!(V_rcv_buf, V_snd_buf, graph)
    @async begin
        fetch(t_V)
        map(store_recv_data!, V, hold_data_size, V_rcv_buf)
        map(split_and_compress!, partition(A), V, change_sparse, perm_sparse)
        A
    end
end

function assemble_matrix_with_compressed_snd_and_with_int_vector_cache!(I, J, V, rows, cols, f)
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
            left_ptr += 1
            right_ptr -= 1
        end
        Tkey = eltype(key)
        change = collect(Tkey, 1:right_ptr)
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
            change[left_ptr] = right_ptr
            left_ptr += 1
            right_ptr -= 1
        end
        right_ptr, change
    end
    function partition_and_setup_cache_snd!(I, J, V, perm, parts_snd, rows_sa, cols)
        n_hold_data, change_index = quick_sort_partition!(rows_sa, I, J, V)
        n_raw_snd_data = length(I) - n_hold_data
        resize!(perm, n_raw_snd_data)
        sizehint!(perm, n_raw_snd_data)
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
        J_snd_data = similar(I, n_snd_data)
        V_snd_data = similar(V, n_snd_data)
        owner_to_p = Dict(owner => i for (i, owner) in enumerate(parts_snd))
        ghost_to_p = [owner_to_p[owner] for owner in ghost_to_owner(rows_sa)]
        ptrs = zeros(Int32, length(parts_snd) + 1)
        for (_, i, _) in nziterator(compressed_snd_data)
            ptrs[ghost_to_p[i]+1] += 1
        end
        length_to_ptrs!(ptrs)
        for (n, (j, i, v)) in enumerate(nziterator(compressed_snd_data))
            p = ghost_to_p[i]
            index = ptrs[p]
            J_snd_data[index] = j
            V_snd_data[index] = v
            I_snd_data[n] = index
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        precompute_nzindex!(perm, compressed_snd_data, J_raw_snd_data, I_raw_snd_data)
        for (i, p) in enumerate(perm)
            perm[i] = I_snd_data[p]
        end
        ghost_to_global_row = ghost_to_global(rows_sa)
        for (_, i, _) in nziterator(compressed_snd_data)
            p = ghost_to_p[i]
            I_snd_data[ptrs[p]] = ghost_to_global_row[i]
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data, change_index, perm
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
    function split_and_compress!(I, J, V, perm, rows_fa, cols_fa)
        n_own_data, change_index = quick_sort_partition!(cols_fa, J, I, V)
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
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        precompute_nzindex!(perm_own, own_own, I_own_own, J_own_own)
        precompute_nzindex!(perm_ghost, own_ghost, I_own_ghost, J_own_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm), change_index, perm
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size, change_snd, perm_snd = map(partition_and_setup_cache_snd!, I, J, V, I_owner, parts_snd, rows_sa, cols) |> tuple_of_arrays
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
        vals_fa, change_sparse, perm_sparse = map(split_and_compress!, I, J, V, J_owner, rows_fa, cols_fa) |> tuple_of_arrays
        cache = (; graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled), cache
    end
end

function assemble_matrix_with_compressed_snd_and_with_int_vector_cache!(A, V, cache)
    function perm_partition!(V, perm::Vector{T}) where {T}
        for (i, p) in enumerate(perm)
            V[i], V[p] = V[p], V[i]
        end
    end
    function partition_and_setup_cache_snd!(V_snd, V, n_hold_data, change_index, perm)
        perm_partition!(V, change_index)
        snd_index = (n_hold_data+1):lastindex(V)
        V_raw_snd_data = view(V, snd_index)
        V_snd_data = V_snd.data
        fill!(V_snd_data, 0)
        for (p, v) in zip(perm, V_raw_snd_data)
            V_snd_data[p] += v
        end
    end
    function store_recv_data!(V, n_hold_data, V_rcv)
        n_data = n_hold_data + length(V_rcv.data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(A, V, change_index, perm)
        perm_partition!(V, change_index)
        n_own_data = length(change_index)
        is_own = firstindex(V):n_own_data
        is_ghost = (n_own_data+1):lastindex(V)
        V_own_own = view(V, is_own)
        V_own_ghost = view(V, is_ghost)
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        sparse_matrix!(A.blocks.own_own, V_own_own, perm_own; reset=false)
        sparse_matrix!(A.blocks.own_ghost, V_own_ghost, perm_ghost; reset=false)
        return
    end
    graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse = cache
    map(partition_and_setup_cache_snd!, V_snd_buf, V, hold_data_size, change_snd, perm_snd)
    t_V = exchange!(V_rcv_buf, V_snd_buf, graph)
    @async begin
        fetch(t_V)
        map(store_recv_data!, V, hold_data_size, V_rcv_buf)
        map(split_and_compress!, partition(A), V, change_sparse, perm_sparse)
        A
    end
end

function assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache!(I, J, V, rows, cols, f)
    function quick_sort_partition!(part, key, values::Vararg{Any,N}) where {N}
        global_to_own_part = global_to_own(part)
        n_change = 0
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
            n_change += 1
            left_ptr += 1
            right_ptr -= 1
        end
        Tkey = eltype(key)
        change = Vector{Tuple{Tkey,Tkey}}(undef, n_change)
        ptr = firstindex(change)
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
            change[ptr] = (left_ptr, right_ptr)
            ptr += 1
            left_ptr += 1
            right_ptr -= 1
        end
        right_ptr, change
    end
    function partition_and_setup_cache_snd!(I, J, V, perm, parts_snd, rows_sa, cols)
        n_hold_data, change_index = quick_sort_partition!(rows_sa, I, J, V)
        n_raw_snd_data = length(I) - n_hold_data
        resize!(perm, n_raw_snd_data)
        sizehint!(perm, n_raw_snd_data)
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
        for (n, (j, i, v)) in enumerate(nziterator(compressed_snd_data))
            p = ghost_to_p[i]
            index = ptrs[p]
            I_snd_data[n] = index
            J_snd_data[index] = j
            V_snd_data[index] = v
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        precompute_nzindex!(perm, compressed_snd_data, J_raw_snd_data, I_raw_snd_data)
        for (i, p) in enumerate(perm)
            perm[i] = I_snd_data[p]
        end
        ghost_to_global_row = ghost_to_global(rows_sa)
        for (_, i, _) in nziterator(compressed_snd_data)
            p = ghost_to_p[i]
            I_snd_data[ptrs[p]] = ghost_to_global_row[i]
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data, change_index, perm
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
    function split_and_compress!(I, J, V, perm, rows_fa, cols_fa)
        n_own_data, change_index = quick_sort_partition!(cols_fa, J, I, V)
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
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        precompute_nzindex!(perm_own, own_own, I_own_own, J_own_own)
        precompute_nzindex!(perm_ghost, own_ghost, I_own_ghost, J_own_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm), change_index, perm
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size, change_snd, perm_snd = map(partition_and_setup_cache_snd!, I, J, V, I_owner, parts_snd, rows_sa, cols) |> tuple_of_arrays
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
        vals_fa, change_sparse, perm_sparse = map(split_and_compress!, I, J, V, J_owner, rows_fa, cols_fa) |> tuple_of_arrays
        cache = (; graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled), cache
    end
end

function assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache!(A, V, cache)
    function perm_partition!(V, perm::Vector{Tuple{T,T}}) where {T}
        for (i, j) in perm
            V[i], V[j] = V[j], V[i]
        end
    end
    function partition_and_setup_cache_snd!(V_snd, V, n_hold_data, change_index, perm)
        perm_partition!(V, change_index)
        snd_index = (n_hold_data+1):lastindex(V)
        V_raw_snd_data = view(V, snd_index)
        V_snd_data = V_snd.data
        fill!(V_snd_data, 0)
        for (p, v) in zip(perm, V_raw_snd_data)
            V_snd_data[p] += v
        end
    end
    function store_recv_data!(V, n_hold_data, V_rcv)
        n_data = n_hold_data + length(V_rcv.data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(A, V, change_index, perm)
        perm_partition!(V, change_index)
        n_own_data = isempty(change_index) ? 1 : last(change_index)[1]
        is_own = firstindex(V):n_own_data
        is_ghost = (n_own_data+1):lastindex(V)
        V_own_own = view(V, is_own)
        V_own_ghost = view(V, is_ghost)
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        sparse_matrix!(A.blocks.own_own, V_own_own, perm_own; reset=false)
        sparse_matrix!(A.blocks.own_ghost, V_own_ghost, perm_ghost; reset=false)
        return
    end
    graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse = cache
    map(partition_and_setup_cache_snd!, V_snd_buf, V, hold_data_size, change_snd, perm_snd)
    t_V = exchange!(V_rcv_buf, V_snd_buf, graph)
    @async begin
        fetch(t_V)
        map(store_recv_data!, V, hold_data_size, V_rcv_buf)
        map(split_and_compress!, partition(A), V, change_sparse, perm_sparse)
        A
    end
end

function assemble_matrix_with_compressed_snd_and_with_auto_cache!(I, J, V, rows, cols, f)
    function quick_sort_partition!(part, key, values::Vararg{Any,N}) where {N}
        global_to_own_part = global_to_own(part)
        n_change = 0
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
            n_change += 1
            left_ptr += 1
            right_ptr -= 1
        end
        Tkey = eltype(key)
        if right_ptr - n_change - n_change <= 0
            change = collect(Tkey, 1:right_ptr)
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
                change[left_ptr] = right_ptr
                left_ptr += 1
                right_ptr -= 1
            end
        else
            change = Vector{Tuple{Tkey,Tkey}}(undef, n_change)
            ptr = firstindex(change)
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
                change[ptr] = (left_ptr, right_ptr)
                ptr += 1
                left_ptr += 1
                right_ptr -= 1
            end
        end
        right_ptr, change
    end
    function partition_and_setup_cache_snd!(I, J, V, perm, parts_snd, rows_sa, cols)
        n_hold_data, change_index = quick_sort_partition!(rows_sa, I, J, V)
        n_raw_snd_data = length(I) - n_hold_data
        resize!(perm, n_raw_snd_data)
        sizehint!(perm, n_raw_snd_data)
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
        for (n, (j, i, v)) in enumerate(nziterator(compressed_snd_data))
            p = ghost_to_p[i]
            index = ptrs[p]
            I_snd_data[n] = index
            J_snd_data[index] = j
            V_snd_data[index] = v
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        precompute_nzindex!(perm, compressed_snd_data, J_raw_snd_data, I_raw_snd_data)
        for (i, p) in enumerate(perm)
            perm[i] = I_snd_data[p]
        end
        ghost_to_global_row = ghost_to_global(rows_sa)
        for (_, i, _) in nziterator(compressed_snd_data)
            p = ghost_to_p[i]
            I_snd_data[ptrs[p]] = ghost_to_global_row[i]
            ptrs[p] += 1
        end
        rewind_ptrs!(ptrs)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, n_hold_data, change_index, perm
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
    function split_and_compress!(I, J, V, perm, rows_fa, cols_fa)
        n_own_data, change_index = quick_sort_partition!(cols_fa, J, I, V)
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
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        precompute_nzindex!(perm_own, own_own, I_own_own, J_own_own)
        precompute_nzindex!(perm_ghost, own_ghost, I_own_ghost, J_own_ghost)
        rows_perm = local_permutation(rows_fa)
        cols_perm = local_permutation(cols_fa)
        split_matrix(blocks, rows_perm, cols_perm), change_index, perm
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd_buf, J_snd_buf, V_snd_buf, hold_data_size, change_snd, perm_snd = map(partition_and_setup_cache_snd!, I, J, V, I_owner, parts_snd, rows_sa, cols) |> tuple_of_arrays
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
        vals_fa, change_sparse, perm_sparse = map(split_and_compress!, I, J, V, J_owner, rows_fa, cols_fa) |> tuple_of_arrays
        cache = (; graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse)
        assembled = true
        PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled), cache
    end
end

function assemble_matrix_with_compressed_snd_and_with_auto_cache!(A, V, cache)
    function partition_and_setup_cache_snd!(V_snd, V, n_hold_data, change_index, perm)
        perm_partition!(v, change::Vector{T}) where {T} = begin
            for (i, p) in enumerate(change)
                v[i], v[p] = v[p], v[i]
            end
        end
        perm_partition!(v, change::Vector{Tuple{T,T}}) where {T} = begin
            for (i, j) in change
                v[i], v[j] = v[j], v[i]
            end
        end
        perm_partition!(V, change_index)
        snd_index = (n_hold_data+1):lastindex(V)
        V_raw_snd_data = view(V, snd_index)
        V_snd_data = V_snd.data
        fill!(V_snd_data, 0)
        for (p, v) in zip(perm, V_raw_snd_data)
            V_snd_data[p] += v
        end
    end
    function store_recv_data!(V, n_hold_data, V_rcv)
        n_data = n_hold_data + length(V_rcv.data)
        resize!(V, n_data)
        sizehint!(V, n_data)
        rcv_index = (n_hold_data+1):n_data
        V[rcv_index] = V_rcv.data
        return
    end
    function split_and_compress!(A, V, change_index, perm)
        perm_partition!(v, change::Vector{T}) where {T} = begin
            for (i, p) in enumerate(change)
                v[i], v[p] = v[p], v[i]
            end
        end
        perm_partition!(v, change::Vector{Tuple{T,T}}) where {T} = begin
            for (i, j) in change
                v[i], v[j] = v[j], v[i]
            end
        end
        get_n_own_data(change_index::Vector{T}) where {T} = length(change_index)
        get_n_own_data(change_index::Vector{Tuple{T,T}}) where {T} = isempty(change_index) ? 1 : last(change_index)[1]
        perm_partition!(V, change_index)
        n_own_data = get_n_own_data(change_index)
        is_own = firstindex(V):n_own_data
        is_ghost = (n_own_data+1):lastindex(V)
        V_own_own = view(V, is_own)
        V_own_ghost = view(V, is_ghost)
        perm_own = view(perm, is_own)
        perm_ghost = view(perm, is_ghost)
        sparse_matrix!(A.blocks.own_own, V_own_own, perm_own; reset=false)
        sparse_matrix!(A.blocks.own_ghost, V_own_ghost, perm_ghost; reset=false)
        return
    end
    graph, V_snd_buf, V_rcv_buf, hold_data_size, change_snd, perm_snd, change_sparse, perm_sparse = cache
    map(partition_and_setup_cache_snd!, V_snd_buf, V, hold_data_size, change_snd, perm_snd)
    t_V = exchange!(V_rcv_buf, V_snd_buf, graph)
    @async begin
        fetch(t_V)
        map(store_recv_data!, V, hold_data_size, V_rcv_buf)
        map(split_and_compress!, partition(A), V, change_sparse, perm_sparse)
        A
    end
end