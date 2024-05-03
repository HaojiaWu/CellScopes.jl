function sp_dim_plot(sp::VisiumHDObject, anno; x_col::String = "x", y_col::String = "y",
    cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing, img_res = "low",
    anno_color::Union{Nothing, Dict} = nothing, adjust_contrast =1.0, adjust_brightness=0.0,
    pt_bg_color = "gray90", x_lims=nothing, y_lims=nothing,width = 500, height = 500, alpha=1,
    stroke_width=0, stroke_color="black", cell_order::Union{Vector{String}, Nothing}=nothing,
    legend_fontsize = 30, do_legend=false, legend_size = 30 , bg_color = "white"
    )
    anno_df = deepcopy(sp.spmetaData)
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
    scale_factor = get_vs_sf(sp; img_res = img_res)
    anno_df[!, x_col] =  anno_df[!, x_col] .* scale_factor
    anno_df[!, y_col] =  anno_df[!, y_col] .* scale_factor
    if isa(cell_highlight, String)
        cell_highlight = [cell_highlight]
    end
    poly = sp.polygonData .* scale_factor
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
    end
    

    anno_df = DataFrames.transform(anno_df, anno => ByRow(x -> anno_color[x]) => :new_color)
    anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
    plt_color = anno_df.new_color
    if img_res == "high"
        img = deepcopy(sp.imageData.highresImage)
    elseif img_res == "low"
        img = deepcopy(sp.imageData.lowresImage)
    else
        img = deepcopy(sp.imageData.fullresImage)
    end
    max_w = minimum([size(img)[1], Int(round(maximum(anno_df[!, x_col])))])
    max_h = minimum([size(img)[2], Int(round(maximum(anno_df[!, y_col])))])
    if isa(x_lims, Nothing)
        x_lims=[1,max_w]
    end
    if isa(y_lims, Nothing)
        y_lims=[1,max_h]
    end
    img = img[x_lims[1]:x_lims[2], y_lims[1]:y_lims[2]]
    img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
    select_fov = filter([:x, :y] => (x, y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], anno_df)
    select_fov = filter(anno => x -> x ∈ (Set(cell_highlight)), select_fov)
    polygon_num = select_fov.ID
    poly = poly[polygon_num]
    poly = [m .- [x_lims[1]-1 y_lims[1]-1] for m in poly]
    plt_color = plt_color[polygon_num]
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
            for (cell1, color1) in zip(cells, colors)
                select_fov1 = filter(anno => ==(cell1), select_fov)
                MK.scatter!(ax1,[0, 1], color = color1, marker=:rect,
                                strokewidth = 0.5,strokecolor=stroke_color, markersize = legend_size, label = cell1)
            end
            MK.Legend(fig[1, 2], ax1, framecolor=:white, labelsize=legend_fontsize)
        end
        MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
    MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
    MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
    return MK.current_figure()
end

