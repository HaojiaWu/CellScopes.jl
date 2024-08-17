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

function transform_coord(df, t_mat; x_old = :x, y_old = :y, x_new = :new_x, y_new = :new_y)
    coordinates = Matrix(df[:, [x_col, y_col]])
    ones_column = ones(size(df, 1), 1)
    homo_coord = hcat(coordinates, ones_column)
    trans_coord = homo_coord * transpose(t_mat)
    df[!, x_new] = trans_coord[:, 1]
    df[!, y_new] = trans_coord[:, 2]
    return df
end

function generate_hd_segcount(xn_dir, vs_dir; t_mat = nothing, img_lims=nothing)
    cell_seg = read_parquet(xn_dir * "cell_boundaries.parquet")
    umap = CSV.read(xn_dir * "analysis/umap/gene_expression_2_components/projection.csv", DataFrame)
    cell_seg = filter(:cell_id=> ∈(Set(umap.Barcode)), cell_seg)
    cell_seg.x = cell_seg.vertex_x ./ 0.2125
    cell_seg.y = cell_seg.vertex_y ./ 0.2125
    gdf1 = DataFrames.combine(groupby(cell_seg, :cell_id)) do df
        DataFrame(geometry = points_to_polygon(df))
    end
    vs_spot = read_parquet(vs_dir * "binned_outputs/square_002um/spatial/tissue_positions.parquet")
    vs_spot = filter(:in_tissue => !=(0), vs_spot)
    if isa(img_lims, Nothing)
        img_lims = [maximum(vs_spot.pxl_row_in_fullres), maximum(vs_spot.pxl_col_in_fullres)]
    end
    vs_spot = vs_spot[
        (vs_spot.pxl_row_in_fullres .> 0) .& 
        (vs_spot.pxl_row_in_fullres .< img_lims[1]) .& 
        (vs_spot.pxl_col_in_fullres .> 0) .& 
        (vs_spot.pxl_col_in_fullres .< img_lims[2]), :]
    if !isa(t_mat, Nothing)
        vs_spot = transform_coord(vs_spot, t_mat; x_old = :pxl_col_in_fullres, y_old = :pxl_row_in_fullres, x_new=:x, y_new = :y)
    else
        vs_spot.x, vs_spot.y = vs_spot.pxl_col_in_fullres, vs_spot.pxl_row_in_fullres
    end
    points = [(x, y) for (x, y) in zip(vs_spot.x, vs_spot.y)]
    points_df = DataFrame(geometry = points, barcode = vs_spot.barcode, x = vs_spot.x, y=vs_spot.y)
    joined_df = FlexiJoins.innerjoin(
            (points_df, gdf1),
            by_pred(:geometry, GO.within, :geometry)
        )
    molecule = parse_molecule(vs_dir)
    molecule = filter(:barcode => ∈(Set(joined_df.barcode)),  molecule)
    molecule = DataFrames.leftjoin(molecule, joined_df, on = :barcode)
    vs_mtx, gene_names, cell_names = pivot_count(molecule)
    return vs_mtx, gene_names, cell_names
end

function reformat_polygons(xn_dir, t_mat)
    cell_seg = read_parquet(xn_dir * "cell_boundaries.parquet")
    umap = CSV.read(xn_dir * "analysis/umap/gene_expression_2_components/projection.csv", DataFrame)
    cell_seg = filter(:cell_id=> ∈(Set(umap.Barcode)), cell_seg)
    cell_seg.vertex_x = cell_seg.vertex_x ./ 0.2125
    cell_seg.vertex_y = cell_seg.vertex_y ./ 0.2125
    inv_vs_mat = inv(vs_mat)
    cell_seg = transform_coord(cell_seg, inv_vs_mat; x_old = :vertex_x, y_old = :vertex_y, x_new=:vertex_y, y_new = :vertex_x)
    grouped = groupby(cell_seg, :cell_id)
    cell_ids = unique(cell_seg.cell_id)
    poly = Vector{Matrix{Float64}}(undef, length(cell_ids))
    n = length(cell_ids)
    println("Reformatting cell polygons...")
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for idx in 1:length(cell_ids)
        cell_data = grouped[idx]
        cell_1 = Matrix(cell_data[!, 2:end])
        poly[idx] = cell_1
        next!(p)
    end
    println("Cell polygons reformatted!")
    return poly
end
