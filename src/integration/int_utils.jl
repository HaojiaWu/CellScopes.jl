function merge_matrices(matrices::Vector{SparseMatrixCSC{T, Int64}}, row_indices::Vector{Vector{String}}) where T
    all_row_indices = sort(unique(vcat(row_indices...)))
    row_dicts = [Dict(row_indices[i][j] => j for j in 1:length(row_indices[i])) for i in 1:length(row_indices)]
    I, J, V = Int[], Int[], Float64[]
    col_offset = 0
    for (mat_idx, (mat, row_dict)) in enumerate(zip(matrices, row_dicts))
        for (i, row_idx) in enumerate(all_row_indices)
            if haskey(row_dict, row_idx)
                row = row_dict[row_idx]
                for (j, val) in enumerate(mat[row, :])
                    if val != 0
                        push!(I, i)
                        push!(J, col_offset + j)
                        push!(V, val)
                    end
                end
            end
        end
        col_offset += size(mat, 2)
    end
    merged_matrix = sparse(I, J, V, length(all_row_indices), col_offset)
    return merged_matrix, all_row_indices
end

@generated function fill_nothing(::Type{T}) where T
    fields = fieldnames(T)
    Expr(:new, T, (:(nothing) for _ in fields)...)
end