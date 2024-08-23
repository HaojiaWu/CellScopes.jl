function sp_dim_plot(sp::PairedObject; 
        anno::Union{Symbol, String}="cluster", 
        xn_anno::Union{Symbol, String}="cluster", 
        vs_anno::Union{Symbol, String}="cluster", 
        data_use::String = "individual", img_use::String = "vs_img",
        anno_color::Union{Nothing, Dict} = nothing, 
        cell_order::Union{Vector{String}, Nothing}=nothing, 
        xn_cell_order::Union{Vector{String}, Nothing}=nothing, 
        vs_cell_order::Union{Vector{String}, Nothing}=nothing,
        cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing, 
        xn_cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing,
        vs_cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing, 
        pt_bg_color="white",
        stroke_width=0.5,
        stroke_color=:transparent, 
        label_size=50, 
        label_color="black", 
        label_offset=(0,0), 
        do_label=false, 
        alpha::Real = 1, 
        legend_ncol = 1,
        xn_legend_ncol = 1, 
        vs_legend_ncol = 1, 
        plot_img = true,
        x_col = "x", 
        y_col = "y", 
        x_lims = nothing, 
        y_lims = nothing,
        marker_size = 2, 
        bg_color = :white,  
        adjust_contrast = 1.0,
        adjust_brightness = 0.0, 
        legend_size = 30, 
        legend_fontsize = 20, 
        do_legend = false,
        height = 500, 
        width = 500
    )
    if isa(x_lims, Nothing)
        x_lims=(minimum(anno_df[!,x_col])-0.05*maximum(anno_df[!,x_col]),1.05*maximum(anno_df[!,x_col]))
    else 
        x_lims = x_lims
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(anno_df[!,y_col])-0.05*maximum(anno_df[!,y_col]),1.05*maximum(anno_df[!,y_col]))
    else 
        y_lims = y_lims
    end
    if data_use == "cellseg"
        anno_df=deepcopy(sp.spmetaData.cell)
        anno_df[!, anno] = string.(anno_df[!, anno])
        if isa(anno, String)
            anno=Symbol(anno)
        end
        if isa(anno_color, Nothing)
            cell_anno=unique(anno_df[!,anno])
            c_map=Colors.distinguishable_colors(length(cell_anno), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
            c_map = "#" .* hex.(c_map)
            anno_color1=Dict(cell_anno .=> c_map)

        else
            anno_color1 = anno_color
        end
        anno_df=DataFrames.transform(anno_df, anno => ByRow(x -> anno_color1[x]) => :new_color)
        anno_df = filter([x_col, y_col] => (x,y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], anno_df)
        anno_df[!, x_col] = anno_df[!, x_col] .- x_lims[1]
        anno_df[!, y_col] = anno_df[!, y_col] .- y_lims[1]
        anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
        if isa(cell_highlight, String)
            cell_highlight = [cell_highlight]
        end
        all_celltypes = unique(anno_df[!,vs_anno])
        if isa(cell_highlight, Nothing)
            cell_highlight = all_celltypes
        end
        anno_df = filter(anno => ∈(Set(cell_highlight)), anno_df)
        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false);
        ax2 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false)
        if img_use in ["xn_img", "vs_img"]
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
            if plot_img
                MK.image!(ax1, img2)
            end
        end
        if isa(cell_order, Nothing)
            cell_anno=unique(anno_df[!,anno])
        else
            cell_anno = cell_order
        end
        for i in cell_anno
            anno_df3=filter(anno => ==(i), anno_df)
            x_ax = anno_df3[!, x_col]
            y_ax = anno_df3[!, y_col]
            colors = unique(anno_df3.new_color)
            if do_legend
                MK.scatter!(ax2, x_ax , y_ax; strokecolor=stroke_color, visible=false,
                    color=colors[1], strokewidth=0, markersize=legend_size, label=i)
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=colors[1], strokewidth=0, markersize=marker_size, label=i)
            else
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=colors[1], strokewidth=0, markersize=marker_size)
            end
        end

        if do_legend
            MK.Legend(fig[1, 2], ax2, framecolor=:white, labelsize=legend_fontsize, nbanks=legend_ncol)
        end
        if do_label
            for i in cell_anno
                anno_df3 = filter(anno => ==(i), anno_df)
                x_ax = anno_df3[!, x_col]
                y_ax = anno_df3[!, y_col]
                MK.text!(ax1, i, position = (mean(x_ax) - label_offset[1], mean(y_ax) - label_offset[2]),align = (:center, :center),font = "Noto Sans Regular",fontsize = label_size,color = label_color)
            end
        end
        return MK.current_figure()
        
    elseif data_use == "individual"
        ### xenium processing
        anno_df=deepcopy(sp.pairedData.xnObj.spmetaData.cell)
        anno_df[!, xn_anno] = string.(anno_df[!, xn_anno])
        if isa(xn_anno, String)
            anno=Symbol(xn_anno)
        end
        if isa(anno_color, Nothing)
            cell_anno=unique(anno_df[!,xn_anno])
            c_map=Colors.distinguishable_colors(length(cell_anno), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
            c_map = "#" .* hex.(c_map)
            anno_color1=Dict(cell_anno .=> c_map)
        else
            anno_color1 = anno_color
        end
        anno_df=DataFrames.transform(anno_df, xn_anno => ByRow(x -> anno_color1[x]) => :new_color)
        anno_df = filter([x_col, y_col] => (x,y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], anno_df)
        anno_df[!, x_col] = anno_df[!, x_col] .- x_lims[1]
        anno_df[!, y_col] = anno_df[!, y_col] .- y_lims[1]
        anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
        img = deepcopy(sp.pairedData.xnObj.imageData)
        if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
            img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
        end
        img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        if isa(xn_cell_highlight, String)
            xn_cell_highlight = [xn_cell_highlight]
        end
        all_celltypes = unique(anno_df[!,vs_anno])
        if isa(xn_cell_highlight, Nothing)
            xn_cell_highlight = all_celltypes
        end
        anno_df = filter(xn_anno => ∈(Set(xn_cell_highlight)), anno_df)
        # visium processing
        anno_df2 = deepcopy(sp.pairedData.vsObj.spmetaData)
        x_col = Symbol(x_col)
        y_col = Symbol(y_col)
        vs_anno = Symbol(vs_anno)
        rename!(anno_df2, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])
        if isa(anno_df2[!, x_col], Vector{String})
            anno_df2[!, x_col] = Float64.(anno_df2[!, x_col])
        end
        if isa(anno_df2[!, y_col], Vector{String})
            anno_df2[!, y_col] = Float64.(anno_df2[!, y_col])
        end
        if isa(vs_cell_highlight, String)
            vs_cell_highlight = [vs_cell_highlight]
        end
        poly = deepcopy(sp.pairedData.vsObj.polygonData)
        all_celltypes = unique(anno_df2[!,vs_anno])
        if isa(vs_cell_highlight, Nothing)
            vs_cell_highlight = all_celltypes
        end
        other_cells = setdiff(all_celltypes, vs_cell_highlight)
        other_color = Dict(other_cells .=> repeat([pt_bg_color], length(other_cells)))
        if isa(anno_color, Nothing)
            c_map=Colors.distinguishable_colors(length(vs_cell_highlight), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
            c_map = "#" .* hex.(c_map)
            cell_color=Dict(vs_cell_highlight .=> c_map)
            anno_color = merge(cell_color, other_color)
        else
            cell_color = anno_color
            anno_color = merge(anno_color, other_color)
        end
        anno_df2 = DataFrames.transform(anno_df2, vs_anno => ByRow(x -> anno_color[x]) => :new_color)
        anno_df2.new_color = [(i, alpha) for i in anno_df2.new_color]
        plt_color = anno_df2.new_color
        img_vs = deepcopy(sp.pairedData.vsObj.imageData.fullresImage)
        min_w = maximum([1, Int(round(minimum(anno_df2[!, x_col])))])
        min_h = maximum([1, Int(round(minimum(anno_df2[!, y_col])))])    
        max_w = minimum([size(img)[1], Int(round(maximum(anno_df2[!, x_col])))])
        max_h = minimum([size(img)[2], Int(round(maximum(anno_df2[!, y_col])))])
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
        img_vs = img_vs[x_lims[1]:x_lims[2], y_lims[1]:y_lims[2]]
        img_vs2 = augment(img_vs, ColorJitter(adjust_contrast, adjust_brightness))

        select_fov = filter([:x, :y] => (x, y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], anno_df2)
        select_fov = filter(vs_anno => ∈(Set(vs_cell_highlight)), select_fov)
        polygon_num = select_fov.ID
        poly = poly[polygon_num]
        poly = [m .- [x_lims[1]-1 y_lims[1]-1] for m in poly]
        plt_color = select_fov.new_color

        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false);
        ax2 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false)
        if plot_img
            MK.image!(ax1, img2)
        end
        if isa(xn_cell_order, Nothing)
            cell_anno=unique(anno_df[!,xn_anno])
        else
            cell_anno = xn_cell_order
        end
        for i in cell_anno
            anno_df3=filter(xn_anno => ==(i), anno_df)
            x_ax = anno_df3[!, x_col]
            y_ax = anno_df3[!, y_col]
            colors = unique(anno_df3.new_color)
            if do_legend
                MK.scatter!(ax2, x_ax , y_ax; strokecolor=stroke_color, visible=false,
                    color=colors[1], strokewidth=0, markersize=legend_size, label=i)
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=colors[1], strokewidth=0, markersize=marker_size, label=i)
            else
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=colors[1], strokewidth=0, markersize=marker_size)
            end
        end
        if do_legend
            MK.Legend(fig[1, 2], ax2, framecolor=:white, labelsize=legend_fontsize, nbanks=xn_legend_ncol)
            ax3 = MK.Axis(fig[1,3]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
        else
            ax3 = MK.Axis(fig[1,2]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
        end
        if do_label
            for i in cell_anno
                anno_df3 = filter(anno => ==(i), anno_df)
                x_ax = anno_df3[!, x_col]
                y_ax = anno_df3[!, y_col]
                MK.text!(ax1, i, position = (mean(x_ax) - label_offset[1], mean(y_ax) - label_offset[2]),align = (:center, :center),font = "Noto Sans Regular",fontsize = label_size,color = label_color)
            end
        end
        if plot_img
            MK.image!(ax3, img_vs2)
        end
        if do_legend
            if !isa(vs_cell_order, Nothing)
            cells = vs_cell_order
            else
            cells = string.(collect(keys(cell_color)))
            end
            colors = [cell_color[i] for i in cells]
            for (cell1, color1) in zip(cells, colors)
                MK.scatter!(ax3,[NaN], [NaN], color = color1, marker=:rect,
                                strokewidth = 0.5,strokecolor=stroke_color, markersize = legend_size, label = cell1)
            end
            MK.Legend(fig[1, 4], ax3, framecolor=:white, labelsize=legend_fontsize, nbanks=vs_legend_ncol)
        end
        MK.poly!(ax3, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
        MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
        MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
        return MK.current_figure()
    else
        error("data_use can only be cellseg or individual.")
    end
end