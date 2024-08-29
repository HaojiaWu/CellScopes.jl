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
    coordinates = Matrix(df[:, [x_old, y_old]])
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
    cell_count = RawCountObject(vs_mtx, cell_names, gene_names)
    return cell_count, molecule
end

function reformat_polygons(xn_dir, t_mat)
    cell_seg = read_parquet(xn_dir * "cell_boundaries.parquet")
    umap = CSV.read(xn_dir * "analysis/umap/gene_expression_2_components/projection.csv", DataFrame)
    cell_seg = filter(:cell_id=> ∈(Set(umap.Barcode)), cell_seg)
    cell_seg.vertex_x = cell_seg.vertex_x ./ 0.2125
    cell_seg.vertex_y = cell_seg.vertex_y ./ 0.2125
    inv_vs_mat = inv(t_mat)
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
    poly = [mat[:, 1:2] for mat in poly]
    println("Cell polygons reformatted!")
    return poly
end

function normalize_paired_object(sp::PairedObject; kwargs...)
    sp = normalize_object(sp; kwargs...)
    sp = normalize_object(sp.pairedData.xnObj; kwargs...)
    sp = normalize_object(sp.pairedData.vsObj; kwargs...)
    return sp
end

### interpolate the out-of-bound pixels arising from rounding errors with the nearby color
function smoothe_img!(img)
    h, w = size(img)
    for y in 2:h-1, x in 2:w-1
        if img[y, x] == RGB{N0f8}(1.0, 1.0, 1.0)  
            neighbors = [img[y-1, x], img[y+1, x], img[y, x-1], img[y, x+1]]
            valid_neighbors = neighbors[neighbors .!= RGB{N0f8}(1.0, 1.0, 1.0)]
            if length(valid_neighbors) > 0
                img[y, x] = mean(valid_neighbors)
            end
        end
    end
end

function process_xn_transcript_data(xn_obj, gene_list;
    x_lims = nothing,
    y_lims = nothing,
    x_col = "x",  
    y_col = "y",  
    gene_colors = nothing,
    bg_tx = false
)
    if isa(gene_list, String)
        gene_list = [gene_list]
    end
    if isa(gene_colors, String)
        gene_colors = [gene_colors]
    end
    df_plt=deepcopy(xn_obj.spmetaData.molecule)
    df_plt.gene = String.(df_plt.gene)
    if isa(x_lims, Nothing)
        x_lims=(minimum(df_plt[!, x_col])-0.05*maximum(df_plt[!, x_col]),1.05*maximum(df_plt[!, x_col]))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(df_plt[!, y_col])-0.05*maximum(df_plt[!, y_col]),1.05*maximum(df_plt[!, y_col]))
    end
    if isa(gene_colors, Nothing)
        c_map=Colors.distinguishable_colors(length(gene_list), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        gene_colors = "#" .* hex.(c_map)
    end
    gene_color=Dict(gene_list .=> gene_colors)
    gene_color["others"] = "gray95"
    from = collect(keys(gene_color))
    to = collect(values(gene_color))
    gene_set = Set(from)
    df_plt.new_gene = ifelse.(in.(df_plt.gene, Ref(gene_set)), df_plt.gene, "others")
    df_plt = map_values(df_plt, :new_gene, :forcolor, from, to)
    df_plt.new_gene = string.(df_plt.new_gene)
    df_plt.forcolor = [(i, alpha) for i in df_plt.forcolor]
    if bg_tx
        all_genes = ["others"; from[from .!= "others"]]
        all_colors = ["gray95"; to[to .!= "gray95"]]
    else
        all_genes = gene_list
        all_colors = gene_colors
    end
    df_plt = @views filter([x_col, y_col] => (x,y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], df_plt)
    df_plt[!, x_col] = df_plt[!, x_col] .- x_lims[1]
    df_plt[!, y_col] = df_plt[!, y_col] .- y_lims[1]
    return df_plt, all_genes, all_colors
end

function process_hd_dimplot_data(hd_obj;
    anno::Union{Symbol, String}="cluster", 
    anno_color::Union{Vector{String}, Nothing} = nothing,
    x_col = "x",  
    y_col = "y", 
    cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing,
    x_lims = nothing, 
    y_lims = nothing,
    pt_bg_color = "transparent",
    alpha::Real = 0.5,
    adjust_contrast= 1.0,
    adjust_brightness = 0.0
)
    anno_df = deepcopy(hd_obj.spmetaData)
    x_col = Symbol(x_col)
    y_col = Symbol(y_col)
    anno = Symbol(anno)
    rename!(anno_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])
    if isa(anno_df[!, x_col], Vector{String})
        anno_df[!, x_col] = Float64.(anno_df[!, x_col])
    end
    if isa(anno_df[!, y_col], Vector{String})
        anno_df[!, y_col] = Float64.(anno_df[!, y_col])
    end
    if isa(cell_highlight, String)
        cell_highlight = [cell_highlight]
    end
    poly = deepcopy(hd_obj.polygonData)
    all_celltypes = unique(anno_df[!,anno])
    if isa(cell_highlight, Nothing)
        cell_highlight = all_celltypes
    end
    other_cells = setdiff(all_celltypes, cell_highlight)
    other_color = Dict(other_cells .=> repeat([pt_bg_color], length(other_cells)))
    if isa(anno_color, Nothing)
        c_map=Colors.distinguishable_colors(length(cell_highlight), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        c_map = "#" .* hex.(c_map)
        cell_color=Dict(cell_highlight .=> c_map)
        anno_color = merge(cell_color, other_color)
    else
        anno_color = merge(anno_color, other_color)
    end
    anno_df = DataFrames.transform(anno_df, anno => ByRow(x -> anno_color[x]) => :new_color)
    anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
    plt_color = anno_df.new_color
    img = deepcopy(hd_obj.imageData.fullresImage)
    min_w = maximum([1, Int(round(minimum(anno_df[!, x_col])))])
    min_h = maximum([1, Int(round(minimum(anno_df[!, y_col])))])    
    max_w = minimum([size(img)[1], Int(round(maximum(anno_df[!, x_col])))])
    max_h = minimum([size(img)[2], Int(round(maximum(anno_df[!, y_col])))])
    if isa(x_lims, Nothing)
        x_lims=[min_w, max_w]
    else
        xlim1 = Int(round(x_lims[1]))
        xlim2 = Int(round(x_lims[2] ))
        x_lims = [xlim1, xlim2]
    end
    if isa(y_lims, Nothing)
        y_lims=[min_h, max_h]
    else
        ylim1 = Int(round(y_lims[1]))
        ylim2 = Int(round(y_lims[2] ))
        y_lims = [ylim1, ylim2]
    end
    img = img[x_lims[1]:x_lims[2], y_lims[1]:y_lims[2]]
    img = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
    select_fov = filter([:x, :y] => (x, y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], anno_df)
    select_fov = filter(anno => ∈(Set(cell_highlight)), select_fov)
    polygon_num = select_fov.ID
    poly = poly[polygon_num]
    poly = [m .- [x_lims[1]-1 y_lims[1]-1] for m in poly]
    plt_color = select_fov.new_color
    return img, poly, cell_color, plt_color
end

function process_hd_featureplot_data(hd_obj, gene;
    color_keys::Union{Vector{String}, Tuple{String}}=["gray94","lemonchiffon","orange","red3"],
    x_col = "x",  
    y_col = "y", 
    hd_layer = "8_um",
    clip = 0.0,
    scale =  false,
    x_lims = nothing, 
    y_lims = nothing,
    adjust_contrast= 1.0,
    adjust_brightness = 0.0
)

    hd_obj = set_default_layer(hd_obj; layer_slot = hd_layer)
    if isa(hd_obj.normCount, Nothing)
        hd_obj = normalize_object(hd_obj)
    end
    norm_count=hd_obj.normCount
    anno_df = deepcopy(hd_obj.spmetaData)
    x_col = Symbol(x_col)
    y_col = Symbol(y_col)
    rename!(anno_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])

    if isa(anno_df[!, x_col], Vector{String})
        anno_df[!, x_col] = Float64.(anno_df[!, x_col])
    end
    if isa(anno_df[!, y_col], Vector{String})
        anno_df[!, y_col] = Float64.(anno_df[!, y_col])
    end
    poly = deepcopy(hd_obj.polygonData)
    c_map = ColorSchemes.ColorScheme([parse(Colorant, color_keys[1]),parse(Colorant, color_keys[2]),parse(Colorant, color_keys[3]),parse(Colorant, color_keys[4])])
    img = deepcopy(hd_obj.imageData.fullresImage)
    min_w = maximum([1, Int(round(minimum(anno_df[!, x_col])))])
    min_h = maximum([1, Int(round(minimum(anno_df[!, y_col])))])    
    max_w = minimum([size(img)[1], Int(round(maximum(anno_df[!, x_col])))])
    max_h = minimum([size(img)[2], Int(round(maximum(anno_df[!, y_col])))])
    if isa(x_lims, Nothing)
        x_lims=[min_w, max_w]
    else
        xlim1 = Int(round(x_lims[1]))
        xlim2 = Int(round(x_lims[2] ))
        x_lims = [xlim1, xlim2]
    end
    if isa(y_lims, Nothing)
        y_lims=[min_h, max_h]
    else
        ylim1 = Int(round(y_lims[1]))
        ylim2 = Int(round(y_lims[2] ))
        y_lims = [ylim1, ylim2]
    end
    img = img[x_lims[1]:x_lims[2], y_lims[1]:y_lims[2]]
    img = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
    poly = [m .- [x_lims[1]-1 y_lims[1]-1] for m in poly]

    gene_expr = subset_count(norm_count; genes = [gene])
    gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
    anno_df.gene = gene_expr
    select_fov = filter([:x, :y, :gene] => (x, y, gene) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2] && gene > clip, anno_df)
    polygon_num = select_fov.ID
    poly2 = poly[polygon_num]
    gene_expr = select_fov.gene
    if scale
        gene_expr = unit_range_scale(gene_expr)
    end
    colors = get(c_map, gene_expr, :extrema)
    plt_color="#" .* hex.(colors)
    return img, poly2, gene_expr, plt_color, c_map
end

function process_paired_featureplot_data(sp::PairedObject, gene::String;
    color_keys::Union{Vector{String}, Tuple{String}}=["gray94","lemonchiffon","orange","red3"],
    x_col = "x",  
    y_col = "y", 
    clip = 0.0,
    scale =  false,
    x_lims = nothing, 
    y_lims = nothing,
    adjust_contrast= 1.0,
    adjust_brightness = 0.0,
    img_use = "xn_img"
)

if isa(sp.normCount, Nothing)
    sp = normalize_object(sp)
end
norm_count=sp.normCount
anno_df = deepcopy(sp.spmetaData.cell)
x_col = Symbol(x_col)
y_col = Symbol(y_col)
if isa(anno_df[!, x_col], Vector{String})
    anno_df[!, x_col] = Float64.(anno_df[!, x_col])
end
if isa(anno_df[!, y_col], Vector{String})
    anno_df[!, y_col] = Float64.(anno_df[!, y_col])
end
poly = deepcopy(sp.polygonData)
c_map = ColorSchemes.ColorScheme([parse(Colorant, color_keys[1]),parse(Colorant, color_keys[2]),parse(Colorant, color_keys[3]),parse(Colorant, color_keys[4])])
if img_use == "xn_img"
    img = deepcopy(sp.pairedData.xnObj.imageData)
elseif img_use == "vs_img"
    img = deepcopy(sp.pairedData.vsObj.imageData.fullresImage)
else
    error("img_use can only be vs_img or xn_img!")
end
if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
    img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
end
img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
poly = [m .- [x_lims[1]-1 y_lims[1]-1] for m in poly]

gene_expr = subset_count(norm_count; genes = [gene])
gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
anno_df = reorder(anno_df, :cell, norm_count.cell_name)
anno_df.gene = gene_expr
select_fov = filter([:x, :y, :gene] => (x, y, gene) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2] && gene > clip, anno_df)
poly_df = deepcopy(sp.spmetaData.polygon)
cell_use = String.(select_fov.cell)
poly_df = filter(:mapped_cell => ∈(Set(cell_use)), poly_df)
polygon_num = poly_df.polygon_number
poly2 = poly[polygon_num]
gene_expr = select_fov.gene
if scale
    gene_expr = unit_range_scale(gene_expr)
end
colors = get(c_map, gene_expr, :extrema)
plt_color="#" .* hex.(colors)
return img2, poly2, gene_expr, plt_color, c_map
end