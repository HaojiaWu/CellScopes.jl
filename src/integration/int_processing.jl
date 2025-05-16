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
    sp.pairedData.xnObj = normalize_object(sp.pairedData.xnObj; kwargs...)
    sp.pairedData.vsObj = normalize_object(sp.pairedData.vsObj; kwargs...)
    return sp
end

### interpolate the out-of-bound pixels arising from rounding errors with the nearby color
function smoothe_img(img)
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
    return img
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
    all_genes = xn_obj.rawCount.gene_name
    if !all(x -> x in all_genes, gene_list)
        error("Some genes are missing. Please make sure that all genes in the gene_list are present in the Xenium data!")
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
    cell_shape = "point",
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
        cell_color=Dict(cell_highlight .=> anno_color)
        anno_color = merge(cell_color, other_color)
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
    if cell_shape != "point"
        return img, poly, cell_color, plt_color
    else
        return img, select_fov, cell_color, plt_color
    end
end

function process_hd_featureplot_data(hd_obj, gene;
    color_keys::Union{Vector{String}, Tuple{String}}=["gray94","lemonchiffon","orange","red3"],
    x_col = "x",  
    y_col = "y", 
    hd_layer = "8_um",
    clip = 0.0,
    x_lims = nothing, 
    y_lims = nothing,
    cell_shape = "bin",
    scale = false,
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
    if cell_shape != "point"
        return img, poly2, gene_expr, plt_color, c_map
    else
        return img, select_fov, gene_expr, plt_color, c_map
    end
end

function process_paired_featureplot_data(sp::PairedObject, gene::String;
    color_keys::Union{Vector{String}, Tuple{String}}=["gray94","lemonchiffon","orange","red3"],
    x_col = "x",  
    y_col = "y", 
    clip = 0.0,
    x_lims = nothing, 
    y_lims = nothing,
    adjust_contrast= 1.0,
    adjust_brightness = 0.0,
    img_use = "xn_img",
    scale = false,
    cell_shape = "point"
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
    if cell_shape == "point"
        return img2, select_fov, gene_expr, plt_color, c_map
    else
        return img2, poly2, gene_expr, plt_color, c_map
    end
end

function process_xn_dimplot_data(sp;
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
    adjust_brightness = 0.0,
    cell_shape = "point"
)
    if cell_shape == "point"
        anno_df=deepcopy(sp.spmetaData.cell)
        anno_df[!, anno] = string.(anno_df[!, anno])
        if isa(x_lims, Nothing)
            x_lims=(minimum(anno_df[!,x_col])-0.05*maximum(anno_df[!,x_col]),1.05*maximum(anno_df[!,x_col]))
        end
        if isa(y_lims, Nothing)
            y_lims=(minimum(anno_df[!,y_col])-0.05*maximum(anno_df[!,y_col]),1.05*maximum(anno_df[!,y_col]))
        end
        if isa(anno, String)
            anno=Symbol(anno)
        end
    
        if isa(cell_highlight, String)
            cell_highlight = [cell_highlight]
        end
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
            c_map = anno_color
            cell_color=Dict(cell_highlight .=> c_map)
            anno_color = merge(cell_color, other_color)
        end
    
        anno_df=DataFrames.transform(anno_df, anno => ByRow(x -> anno_color[x]) => :new_color)
        anno_df = filter(anno => ∈(Set(cell_highlight)), anno_df)
        anno_df = filter([x_col, y_col] => (x,y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], anno_df)
        anno_df[!, x_col] = anno_df[!, x_col] .- x_lims[1]
        anno_df[!, y_col] = anno_df[!, y_col] .- y_lims[1]
        anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
        if isa(sp, PairedObject)
            img = deepcopy(sp.pairedData.xnObj.imageData)
        else
            img = deepcopy(sp.imageData)
        end
        if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
            img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
        end
        img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        return img2, anno_df, c_map
    elseif cell_shape == "polygon"
        if isa(x_lims, Nothing)
            x_lims=(minimum(sp.spmetaData.cell.x)-0.05*maximum(sp.spmetaData.cell.x),1.05*maximum(sp.spmetaData.cell.x))
        end
        if isa(y_lims, Nothing)
            y_lims=(minimum(sp.spmetaData.cell.y)-0.05*maximum(sp.spmetaData.cell.y),1.05*maximum(sp.spmetaData.cell.y))
        end
        if isa(cell_highlight, String)
            cell_highlight = [cell_highlight]
        end
        if isa(anno_color, String)
            anno_color = [anno_color]
        end
        anno_df=deepcopy(sp.spmetaData.polygon)
        polygons=deepcopy(sp.polygonData)
        if isa(anno, String)
            anno=Symbol(anno)
        end
        all_celltypes = unique(anno_df[!,anno])
        if isa(cell_highlight, Nothing)
            cell_highlight = all_celltypes
        end
        other_cells = setdiff(all_celltypes, cell_highlight)
        other_color = Dict(other_cells .=> repeat([pt_bg_color], length(other_cells)))
        if isa(anno_color, Nothing)
            c_map= Colors.distinguishable_colors(length(cell_highlight), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
            c_map = "#" .* hex.(c_map)
            cell_color=Dict(cell_highlight .=> c_map)
            anno_color = merge(cell_color, other_color)
        else
            if length(cell_highlight) !== length(anno_color)
                error("The number of colors must equal to the number of cell types!")
            end
            c_map = anno_color
            cell_color=Dict(cell_highlight .=> c_map)
            anno_color = merge(cell_color, other_color)
        end
        anno_df = DataFrames.transform(anno_df, anno => ByRow(x -> anno_color[x]) => :new_color)
        anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
        plt_color = anno_df.new_color
        select_fov = filter([:x, :y] => (x, y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], sp.spmetaData.cell)
        subset_poly = filter(:mapped_cell => ∈(Set(select_fov.cell)), sp.spmetaData.polygon)
        subset_poly = filter(anno => ∈(Set(cell_highlight)), subset_poly)
        polygon_num = subset_poly.polygon_number
        polygons = polygons[polygon_num]
        plt_color1 = plt_color[polygon_num]
        polygons = [m .- [x_lims[1]-1 y_lims[1]-1] for m in polygons]
        if isa(sp, PairedObject)
            img = deepcopy(sp.pairedData.xnObj.imageData)
        else
            img = deepcopy(sp.imageData)
        end
        if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
            img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
        end
        img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        return img2, polygons, cell_color, plt_color1, c_map
    else
        error("""cell_shape can only be "point" or "polygon"!""")
    end
end

function get_affine_matrix(degrees::Union{Float64, Int64}, center::Union{Tuple{Float64, Float64}, Vector{Float64}}; clockwise = true)
    θ = deg2rad(degrees)
    c, s = cos(θ), sin(θ)
    if clockwise
        s = -s
    end
    x0, y0 = center
    return c, -s, x0 * (1 - c) + y0 * s,
           s,  c, y0 * (1 - c) - x0 * s
end


function rotate_coord(df::DataFrame, degrees::Union{Float64, Int64}, center::Union{Tuple{Float64, Float64}, Vector{Float64}};
    x_col::Union{Symbol, String} = :x,
    y_col::Union{Symbol, String} = :y,
    clockwise = true)
    x_col = Symbol(x_col)
    y_col = Symbol(y_col)
    x = df[!, x_col]
    y = df[!, y_col]
    a11, a12, a13, a21, a22, a23 = get_affine_matrix(degrees, center; clockwise=clockwise)
    x_rot = Vector{Float64}(undef, length(x))
    y_rot = Vector{Float64}(undef, length(y))
    @. x_rot = a11 * x + a12 * y + a13
    @. y_rot = a21 * x + a22 * y + a23
    return DataFrame(x = x_rot, y = y_rot)
end


function rotate_object(sp::Union{XeniumObject, VisiumHDObject}, degree; 
    center = nothing, clockwise = true, x_col = :x, y_col = :y)

    if isa(sp, XeniumObject)
        img = deepcopy(sp.imageData)
        df_plt = deepcopy(sp.spmetaData.cell)
    elseif isa(sp, VisiumHDObject)
        img = deepcopy(sp.imageData.fullresImage)
        df_plt = deepcopy(sp.spmetaData)
        rename!(df_plt, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, :x, :y])
        df_plt = df_plt[!, [:cell, :x, :y]]
    else
        error("sp can only be a Xenium or VisiumHD object.")
    end

    if center === nothing
        center = [mean(df_plt[!, x_col]), mean(df_plt[!, y_col])]
    end
    df_img, colors = img_to_df(img)
    println("\033[1;34mRotating image...\033[0m") 
    rot_coords = rotate_coord(df_img, degree, center; clockwise = clockwise)
    println("Image was rotated!")
    rot_x = rot_coords[:, 1]
    rot_y = rot_coords[:, 2]
    println("\033[1;34mRotating cell coordinates...\033[0m") 
    rot_df = rotate_coord(df_plt, degree, center; clockwise = clockwise)
    println("Coordinates were rotated!")
    x_offset = minimum(rot_x)
    y_offset = minimum(rot_y)
    rot_x .-= x_offset
    rot_y .-= y_offset
    rot_df.x .-= x_offset
    rot_df.y .-= y_offset
    new_img = df_to_img(rot_x, rot_y, colors)
    println("\033[1;34mSmoothing image...\033[0m") 
    new_img = smoothe_img(new_img)
    println("Image was smoothed!")
    min_h = Int64(round(minimum(rot_df.y) / 2))
    min_w = Int64(round(minimum(rot_df.x) / 2))
    new_img2 = new_img[(min_h+1):(end - min_h), (min_w+1):(end - min_w)]
    println("\033[1;34mRemoving the outlier pixels...\033[0m") 
    flip_bg_color!(new_img2)
    println("All done!")
    rot_df.x .-= min_h
    rot_df.y .-= min_w
    return new_img2, rot_df
end

function rotate_paired_object(sp::PairedObject, degree; 
    center=nothing, clockwise = true, x_col = :x, y_col = :y)
    @info "Graph rotation in progress. This process updates all data within the paired objects, including images, coordinates, counts, and cell boundaries. Please allow some time for completion."
    xn_img = deepcopy(sp.pairedData.xnObj.imageData)
    xn_df = deepcopy(sp.pairedData.xnObj.spmetaData.cell)
    vs_img = deepcopy(sp.pairedData.vsObj.imageData.fullresImage)
    vs_df = deepcopy(sp.pairedData.vsObj.spmetaData)
    rename!(vs_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, :x, :y])
    #vs_df = vs_df[!, [:cell, :x, :y]]
    if center === nothing
        center = [mean(xn_df[!, x_col]), mean(xn_df[!, y_col])]
    end

    xn_coord, colors_xn = img_to_df(xn_img)
    println("\033[1;34mRotating the xenium image...\033[0m") 
    xn_rot_coord = rotate_coord(xn_coord, degree, center; clockwise = clockwise)
    println("Xenium image was rotated!")
    xn_rot_x = xn_rot_coord[:, 1]
    xn_rot_y = xn_rot_coord[:, 2]
    println("\033[1;34mRotating the xenium cell coordinates...\033[0m") 
    xn_rot_df = rotate_coord(xn_df, degree, center; clockwise = clockwise)
    xn_rot_df.cell = xn_df.cell
    println("Xenium coordinates were rotated!")
    x_offset1 = minimum(xn_rot_x)
    y_offset1 = minimum(xn_rot_y)
    vs_coord, colors_vs = img_to_df(vs_img)
    println("\033[1;34mRotating the visium image...\033[0m") 
    vs_rot_coord = rotate_coord(vs_coord, degree, center; clockwise = clockwise)
    println("Visium image was rotated!")
    vs_rot_x = vs_rot_coord[:, 1]
    vs_rot_y = vs_rot_coord[:, 2]
    println("\033[1;34mRotating the visium cell coordinates...\033[0m") 
    vs_rot_df = rotate_coord(vs_df, degree, center; clockwise = clockwise)
    vs_rot_df.cell = vs_df.cell
    println("Visium coordinates were rotated!")
    x_offset2 = minimum(vs_rot_x)
    y_offset2 = minimum(vs_rot_y)
    x_offset = x_offset1 <= x_offset2 ? x_offset1 : x_offset2
    y_offset = y_offset1 <= y_offset2 ? y_offset1 : y_offset2
    xn_rot_x .-= x_offset
    xn_rot_y .-= y_offset
    xn_rot_df.x .-= x_offset
    xn_rot_df.y .-= y_offset
    new_img_xn = df_to_img(xn_rot_x, xn_rot_y, colors_xn)
    println("\033[1;34mSmoothing the xenium image...\033[0m") 
    new_img_xn = smoothe_img(new_img_xn)
    println("Xenium image was smoothed!")
    flip_bg_color!(new_img_xn)
    vs_rot_x .-= x_offset
    vs_rot_y .-= y_offset
    vs_rot_df.x .-= x_offset
    vs_rot_df.y .-= y_offset
    new_img_vs = df_to_img(vs_rot_x, vs_rot_y, colors_vs)
    println("\033[1;34mSmoothing the visium image...\033[0m") 
    new_img_vs = smoothe_img(new_img_vs)
    println("Visium image was smoothed!")
    flip_bg_color!(new_img_vs)
    println("\033[1;34mFinal images and coordinates trimming...\033[0m") 
    min_h = minimum(xn_rot_df.y) - maximum(xn_rot_df.y) * 0.05
    min_w = minimum(xn_rot_df.x) - maximum(xn_rot_df.x) * 0.05
    min_h = Int(round(min_h))
    min_w = Int(round(min_w))
    ref_img = deepcopy(new_img_xn)
    new_img_xn = new_img_xn[min_w:(size(ref_img)[1]-min_w), min_h:(size(ref_img)[2]-min_h)];
    new_img_vs = new_img_vs[min_w:(size(ref_img)[1]-min_w), min_h:(size(ref_img)[2]-min_h)];
    xn_rot_df = filter([:x, :y] => (x,y) -> min_w < x < (size(ref_img)[1]-min_w) && min_h < y < (size(ref_img)[2]-min_h), xn_rot_df)
    vs_rot_df = filter([:x, :y] => (x,y) -> min_w < x < (size(ref_img)[1]-min_w) && min_h < y < (size(ref_img)[2]-min_h), vs_rot_df);
    xn_rot_df.x .-= min_w
    xn_rot_df.y .-= min_h
    vs_rot_df.x .-= min_w
    vs_rot_df.y .-= min_h
    sp.pairedData.xnObj.imageData = new_img_xn
    cell_set = Set(xn_rot_df.cell)
    xn_df = filter(:cell => ∈(cell_set), xn_df)
    xn_df.x = xn_rot_df.x
    xn_df.y = xn_rot_df.y
    sp.pairedData.xnObj.spmetaData.cell = xn_df
    sp.spmetaData.cell = xn_df
    sp.pairedData.vsObj.imageData.fullresImage = new_img_vs
    cell_set = Set(vs_rot_df.cell)
    vs_df = filter(:cell => ∈(cell_set), vs_df)
    vs_df.x = vs_rot_df.x
    vs_df.y = vs_rot_df.y
    rename!(vs_df, [:cell, :x, :y] .=> [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] )
    sp.pairedData.vsObj.spmetaData = vs_df
    if !isdefined(sp, :normCount)
        println("Normalizing data...")
        sp = normalize_paired_object(sp)
    end
    cell_filtered = filter(:cell => ∈(Set(sp.normCount.cell_name)), sp.spmetaData.cell)
    sp.spmetaData.cell = cell_filtered
    sp.pairedData.xnObj.spmetaData.cell = cell_filtered
    sp.pairedData.vsObj.normCount = subset_count(sp.pairedData.vsObj.normCount; cells=sp.pairedData.vsObj.spmetaData.barcode)
    @info "All done!"
    return sp
end