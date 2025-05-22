function sp_dim_plot(sp::VisiumHDObject, anno; x_col::String = "x", y_col::String = "y",cell_shape = "polygon",
    cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing, img_res = "low",
    anno_color::Union{Nothing, Dict} = nothing, adjust_contrast =1.0, adjust_brightness=0.0,marker_size=2,
    pt_bg_color = "gray90", x_lims=nothing, y_lims=nothing,width = 500, height = 500, alpha=1,
    stroke_width=0, stroke_color="black", cell_order::Union{Vector{String}, Nothing}=nothing,
    legend_fontsize = 30, do_legend=false, legend_size = 30 , bg_color = "white"
    )
    if isa(sp.alterImgData, Nothing)
        anno_df = deepcopy(sp.spmetaData)
        x_col = Symbol(x_col)
        y_col = Symbol(y_col)
        anno = Symbol(anno)
        rename!(anno_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])
    else
        anno_df = deepcopy(sp.alterImgData.posData.positions[img_res * "_pos"])
        anno_df[!, anno] = sp.spmetaData[!, anno]
    end
    if isa(anno_df[!, x_col], Vector{String})
        anno_df[!, x_col] = Float64.(anno_df[!, x_col])
    end
    if isa(anno_df[!, y_col], Vector{String})
        anno_df[!, y_col] = Float64.(anno_df[!, y_col])
    end
    scale_factor = get_vs_sf(sp; img_res = img_res)
    anno_df[!, x_col] =  anno_df[!, x_col] .* scale_factor
    anno_df[!, y_col] =  anno_df[!, y_col] .* scale_factor
    if isa(cell_highlight, String)
        cell_highlight = [cell_highlight]
    end
    if cell_shape == "polygon"
        if isa(sp.alterImgData, Nothing)
            poly = deepcopy(sp.polygonData)
            poly = poly .* scale_factor
        else
            poly = deepcopy(sp.alterImgData.polyData.polygons[img_res * "_poly"])
            poly = poly .* scale_factor
        end
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
        cell_color = anno_color
        anno_color = merge(anno_color, other_color)
    end
    anno_df = DataFrames.transform(anno_df, anno => ByRow(x -> anno_color[x]) => :new_color)
    anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
    plt_color = anno_df.new_color
    if isa(sp.alterImgData, Nothing)
        if img_res == "high"
            img = deepcopy(sp.imageData.highresImage)
        elseif img_res == "low"
            img = deepcopy(sp.imageData.lowresImage)
        else
            img = deepcopy(sp.imageData.fullresImage)
        end
    else
        img = deepcopy(sp.alterImgData.imgData.imgs[img_res])
    end
    min_w = maximum([1, Int(round(minimum(anno_df[!, x_col])))])
    min_h = maximum([1, Int(round(minimum(anno_df[!, y_col])))])    
    max_w = minimum([size(img)[1], Int(round(maximum(anno_df[!, x_col])))])
    max_h = minimum([size(img)[2], Int(round(maximum(anno_df[!, y_col])))])
    if isa(x_lims, Nothing)
        x_lims=[min_w, max_w]
    else
        xlim1 = Int(round(x_lims[1] * scale_factor))
        xlim2 = Int(round(x_lims[2] * scale_factor))
        x_lims = [xlim1, xlim2]
    end
    if isa(y_lims, Nothing)
        y_lims=[min_h, max_h]
    else
        ylim1 = Int(round(y_lims[1] * scale_factor))
        ylim2 = Int(round(y_lims[2] * scale_factor))
        y_lims = [ylim1, ylim2]
    end
    img = img[x_lims[1]:x_lims[2], y_lims[1]:y_lims[2]]
    img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
    select_fov = deepcopy(anno_df)
    filter!([:x, :y] => (x, y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], select_fov)
    filter!(anno => x -> x ∈ (Set(cell_highlight)), select_fov)
    if cell_shape == "polygon"
        polygon_num = select_fov.ID
        poly = poly[polygon_num]
        poly = [m .- [x_lims[1]-1 y_lims[1]-1] for m in poly]
    end
    plt_color = select_fov.new_color
    fig = MK.Figure(size=(width, height))
    ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
        MK.image!(ax1, img2)
        if do_legend
            if !isa(cell_order, Nothing)
               cells = cell_order
            else
               cells = string.(collect(keys(cell_color)))
            end
            colors = [cell_color[i] for i in cells]
            if cell_shape == "polygon"
                for (cell1, color1) in zip(cells, colors)
                    MK.scatter!(ax1,[NaN], [NaN], color = color1, marker=:rect,
                                    strokewidth = 0.5,strokecolor=stroke_color, markersize = legend_size, label = cell1)
                end
            else
                for (cell1, color1) in zip(cells, colors)
                    MK.scatter!(ax1,[NaN], [NaN], color = color1, marker=:circle,
                                    strokewidth = 0.5,strokecolor=stroke_color, markersize = legend_size, label = cell1)
                end
            end
            MK.Legend(fig[1, 2], ax1, framecolor=:white, labelsize=legend_fontsize)
        end
        if cell_shape == "polygon"
            MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
        else
            select_fov[!, "x"] = select_fov[!, "x"] .- x_lims[1]
            select_fov[!, "y"] = select_fov[!, "y"] .- y_lims[1]
            MK.scatter!(ax1, select_fov.x , select_fov.y; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width,markersize= marker_size)
        end
    MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
    MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
    return fig
end

function sp_feature_plot(sp::VisiumHDObject, gene_list::Union{String, Vector{String}, Tuple{String}};
    color_keys::Union{Vector{String}, Tuple{String}}=["gray94","lemonchiffon","orange","red3"],
    x_lims=nothing, y_lims=nothing,width=900,height=1000,stroke_width=0,stroke_color="black", img_res = "low",
    titlesize::Int64=24, scale::Bool = false, bg_color = "white", x_col::String = "x", y_col::String = "y",
    adjust_contrast =1.0, adjust_brightness=0.0, alpha::Real=1, clip = 0.0
    )
    if isa(gene_list, String)
        gene_list = [gene_list]
    end
    n_rows = Int(ceil(length(gene_list) / 3))
    if length(gene_list) < 4
        n_cols = length(gene_list)
    else
        n_cols = 3
    end
    norm_count=sp.normCount
    if isa(sp.alterImgData, Nothing)
        anno_df = deepcopy(sp.spmetaData)
        x_col = Symbol(x_col)
        y_col = Symbol(y_col)
        #anno = Symbol(anno)
        rename!(anno_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])
    else
        anno_df = deepcopy(sp.alterImgData.posData.positions[img_res * "_pos"])
        #anno_df[!, anno] = sp.spmetaData[!, anno]
    end
    if isa(anno_df[!, x_col], Vector{String})
        anno_df[!, x_col] = Float64.(anno_df[!, x_col])
    end
    if isa(anno_df[!, y_col], Vector{String})
        anno_df[!, y_col] = Float64.(anno_df[!, y_col])
    end
    scale_factor = get_vs_sf(sp; img_res = img_res)
    anno_df[!, x_col] =  anno_df[!, x_col] .* scale_factor
    anno_df[!, y_col] =  anno_df[!, y_col] .* scale_factor
    if isa(sp.alterImgData, Nothing)
        poly = deepcopy(sp.polygonData)
        poly = poly .* scale_factor
    else
        poly = deepcopy(sp.alterImgData.polyData.polygons[img_res * "_poly"])
        poly = poly .* scale_factor
    end
    c_map = ColorSchemes.ColorScheme([parse(Colorant, color_keys[1]),parse(Colorant, color_keys[2]),parse(Colorant, color_keys[3]),parse(Colorant, color_keys[4])])
    if isa(sp.alterImgData, Nothing)
        if img_res == "high"
            img = deepcopy(sp.imageData.highresImage)
        elseif img_res == "low"
            img = deepcopy(sp.imageData.lowresImage)
        else
            img = deepcopy(sp.imageData.fullresImage)
        end
    else
        img = deepcopy(sp.alterImgData.imgData.imgs[img_res])
    end
    min_w = maximum([1, Int(round(minimum(anno_df[!, x_col])))])
    min_h = maximum([1, Int(round(minimum(anno_df[!, y_col])))])    
    max_w = minimum([size(img)[1], Int(round(maximum(anno_df[!, x_col])))])
    max_h = minimum([size(img)[2], Int(round(maximum(anno_df[!, y_col])))])
    if isa(x_lims, Nothing)
        x_lims=[min_w,max_w]
    else
        xlim1 = Int(round(x_lims[1] * scale_factor))
        xlim2 = Int(round(x_lims[2] * scale_factor))
        x_lims = [xlim1, xlim2]
    end
    if isa(y_lims, Nothing)
        y_lims=[min_h,max_h]
    else
        ylim1 = Int(round(y_lims[1] * scale_factor))
        ylim2 = Int(round(y_lims[2] * scale_factor))
        y_lims = [ylim1, ylim2]
    end
    poly = [m .- [x_lims[1]-1 y_lims[1]-1] for m in poly]
    img = img[x_lims[1]:x_lims[2], y_lims[1]:y_lims[2]]
    img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
    fig = MK.Figure(size = (width * n_cols, height * n_rows))
    for (i, gene) in enumerate(gene_list)
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
        plt_color = [(i, alpha) for i in plt_color]
        c_map2 = [(i, alpha) for i in c_map]
        n_row = Int(ceil(i/3))
        if i < 4
            n_col1 = 2i-1
            n_col2 = 2i
        else
            n_col1 = 2*(i-3*(n_row-1))-1
            n_col2 = 2*(i-3*(n_row-1))
        end
        ax1 = MK.Axis(fig[n_row,n_col1]; backgroundcolor = bg_color, xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
        xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false,yreversed=false, title = gene_list[i], 
        titlesize = titlesize, xlabel = "", ylabel = "", 
        xlabelsize = titlesize -4, ylabelsize = titlesize -4)
        MK.image!(ax1, img2)
        MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in poly2]; strokecolor=stroke_color, 
                    color=plt_color, strokewidth=stroke_width,label="")
        MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
        MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
        MK.Colorbar(fig[n_row,n_col2], label = "", colormap = c_map2, width=10, limits = (0, maximum(gene_expr)))
    end
    return fig
end
