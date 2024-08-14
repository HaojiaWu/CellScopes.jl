function run_harmony(sc_obj::get_object_group("All"), batch::Union{String, Symbol}; kwargs...)        
    @info "The run_harmony function is a Julia implementation of the Harmony for data integration. Please read the original paper for the algorithm details: https://www.nature.com/articles/s41592-019-0619-0. The Julia codes are based on a python implementation of Harmony (harmonypy): https://github.com/slowkow/harmonypy"
    pca_mat = sc_obj.dimReduction.pca.cell_embedding
    metadata = sc_obj.metaData
    if isa(batch, String)
        batch = Symbol(batch)
    end
    ho = HarmonyObject(pca_mat, metadata, batch; kwargs...)
    harmony_matrix = Matrix{Float64}(ho.Z_corr')
    sc_obj.dimReduction.pca.cell_embedding = harmony_matrix
    return sc_obj
end

function points_to_polygon(df)
    poly = GI.Polygon([GI.LinearRing([(x, y) for (x, y) in zip(df.x, df.y)])])
    return poly
end

function pivot_count(molecule)
    gene_to_idx = Dict{String, Int}()
    cell_to_idx = Dict{String, Int}()
    gene_names = String[]
    cell_names = String[]
    next_gene_idx = 0
    next_cell_idx = 0
    I = Vector{Int}(undef, size(molecule, 1))
    J = Vector{Int}(undef, size(molecule, 1))
    V = Vector{Float64}(undef, size(molecule, 1))
    idx = 1
    for row in eachrow(molecule)
        gene = row[:gene]
        cell_id = row[:cell_id]
        count = row[:count]
        if haskey(gene_to_idx, gene)
            gene_idx = gene_to_idx[gene]
        else
            next_gene_idx += 1
            gene_idx = next_gene_idx
            gene_to_idx[gene] = gene_idx
            push!(gene_names, gene)
        end
        if haskey(cell_to_idx, cell_id)
            cell_idx = cell_to_idx[cell_id]
        else
            next_cell_idx += 1
            cell_idx = next_cell_idx
            cell_to_idx[cell_id] = cell_idx
            push!(cell_names, cell_id)
        end
        I[idx] = gene_idx
        J[idx] = cell_idx
        V[idx] = count
        idx += 1
    end
    num_genes = next_gene_idx
    num_cells = next_cell_idx
    sparse_mtx = sparse(I, J, V, num_genes, num_cells)
    return sparse_mtx, gene_names, cell_names
end