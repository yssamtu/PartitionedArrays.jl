function matrix_assembly(f, I, J, V, rows, cols)
    function setup_cache_snd_with_cleanup!(I, J, V, I_owner, parts_snd, rows)
        function filter_out_without_resize!(global_to_own_row::AbstractVector, I::AbstractVector, J::AbstractVector, V::AbstractVector)
            index = firstindex(I)
            for (i, j, v) in zip(I, J, V)
                @inbounds I[index] = i
                @inbounds J[index] = j
                @inbounds V[index] = v
                index = ifelse((global_to_own_row[i] != 0)::Bool, nextind(I, index), index)
            end
            index - 1
        end
        global_to_own_row = global_to_own(rows)
        snd_index = findall(i -> global_to_own_row[i] == 0, I)
        owner_data = view(I_owner, snd_index)
        gen = (owner => i for (i, owner) in enumerate(parts_snd))
        owner_to_p = Dict(gen)
        ptrs = zeros(Int32, length(parts_snd) + 1)
        for id in owner_data
            ptrs[owner_to_p[id]+1] += 1
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end] - 1
        Ti = eltype(I)
        Tv = eltype(V)
        I_snd_data = zeros(Ti, ndata)
        J_snd_data = zeros(Ti, ndata)
        V_snd_data = zeros(Tv, ndata)
        for (i, j, v, owner) in zip(I[snd_index], J[snd_index], V[snd_index], I_owner[snd_index])
            p = ptrs[owner_to_p[owner]]
            I_snd_data[p] = i
            J_snd_data[p] = j
            V_snd_data[p] = v
            ptrs[owner_to_p[owner]] += 1
        end
        rewind_ptrs!(ptrs)
        data_size = filter_out_without_resize!(global_to_own_row, I, J, V)
        I_snd = JaggedArray(I_snd_data, ptrs)
        J_snd = JaggedArray(J_snd_data, ptrs)
        V_snd = JaggedArray(V_snd_data, ptrs)
        I_snd, J_snd, V_snd, data_size
    end
    function store_recv_data!(I, J, V, data_size, I_rcv, J_rcv, V_rcv)
        total_size = data_size + length(I_rcv.data)
        resize!(I, total_size)
        sizehint!(I, total_size)
        resize!(J, total_size)
        sizehint!(J, total_size)
        resize!(V, total_size)
        sizehint!(V, total_size)
        rcv_index = (data_size+1):total_size
        I[rcv_index] .= I_rcv.data
        J[rcv_index] .= J_rcv.data
        V[rcv_index] .= V_rcv.data
        I, J, V
    end
    function split_and_compress(I, J, V, rows_fa, cols_fa, f, combine=+)
        global_to_own_col = global_to_own(cols_fa)
        is_ghost = findall(j -> global_to_own_col[j] == 0, J)
        is_own = findall(j -> global_to_own_col[j] != 0, J)
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
        own_own = f(I_own_own, J_own_own, V_own_own, n_own_rows, n_own_cols, combine)
        own_ghost = f(I_own_ghost, J_own_ghost, V_own_ghost, n_own_rows, n_ghost_cols, combine)
        ghost_own = f(Ti[], Ti[], Tv[], n_ghost_rows, n_own_cols, combine)
        ghost_ghost = f(Ti[], Ti[], Tv[], n_ghost_rows, n_ghost_cols, combine)
        blocks = split_matrix_blocks(own_own, own_ghost, ghost_own, ghost_ghost)
        split_matrix(blocks, local_permutation(rows_fa), local_permutation(cols_fa))
    end
    I_owner = find_owner(rows, I)
    rows_sa = map(union_ghost, rows, I, I_owner)
    parts_snd, parts_rcv = assembly_neighbors(rows_sa)
    I_snd, J_snd, V_snd, data_size = map(setup_cache_snd_with_cleanup!, I, J, V, I_owner, parts_snd, rows) |> tuple_of_arrays
    graph = ExchangeGraph(parts_snd, parts_rcv)
    t_I = exchange(I_snd, graph)
    t_J = exchange(J_snd, graph)
    t_V = exchange(V_snd, graph)
    I_rcv = fetch(t_I)
    J_rcv = fetch(t_J)
    V_rcv = fetch(t_V)
    map(store_recv_data!, I, J, V, data_size, I_rcv, J_rcv, V_rcv)
    J_owner = find_owner(cols, J)
    rows_fa = rows
    cols_fa = map(union_ghost, cols, J, J_owner)
    vals_fa = map((I, J, V, rows_fa, cols_fa) -> split_and_compress(I, J, V, rows_fa, cols_fa, f), I, J, V, rows_fa, cols_fa)
    assembled = true
    PSparseMatrix(vals_fa, rows_fa, cols_fa, assembled)
end