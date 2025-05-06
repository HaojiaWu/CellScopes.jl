function sp_dim_plot(sp::PairedObject; 
    anno::Union{Symbol, String}="cluster", 
    xn_anno::Union{Symbol, String}="cluster", 
    vs_anno::Union{Symbol, String}="cluster", 
    data_use::String = "individual", 
    img_use::String = "vs_img",
    anno_color::Union{Vector{String}, Nothing} = nothing, 
    xn_anno_color::Union{Vector{String}, Nothing}= nothing, 
    vs_anno_color::Union{Vector{String}, Nothing} = nothing, 
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
    hd_layer = "8_um",
    marker_size = 2, 
    bg_color = :white,  
    adjust_contrast = 1.0,
    adjust_brightness = 0.0, 
    legend_size = 30, 
    legend_fontsize = 20, 
    do_legend = false,
    height = 500, 
    width = 500,
    xn_cell_shape = "point"
)
    if data_use == "cellseg"
        if xn_cell_shape == "point"
            img, anno_df = process_xn_dimplot_data(sp; anno=anno, anno_color=anno_color, x_col = x_col,  y_col = y_col, 
                cell_highlight= cell_highlight, x_lims = x_lims, y_lims = y_lims, pt_bg_color = pt_bg_color, alpha=alpha,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = xn_cell_shape
            )
        else
            img, polygons, cell_color, plt_color1, c_map = process_xn_dimplot_data(sp; anno= anno, anno_color= anno_color, x_col = x_col,  y_col = y_col, 
                cell_highlight= cell_highlight, x_lims = x_lims, y_lims = y_lims, pt_bg_color = pt_bg_color, alpha=alpha,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = xn_cell_shape
            )
        end
        if img_use == "vs_img"
            img = deepcopy(sp.pairedData.vsObj.imageData.fullresImage)
            if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
                img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
            end
            img = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        end
        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false);
        ax2 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false)
        if plot_img
            MK.image!(ax1, img)
        end
        if xn_cell_shape == "point"
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
            MK.Legend(fig[1, 2], ax1, String(anno), framecolor=:white, labelsize=legend_fontsize, nbanks=legend_ncol, titlesize=20, titlefont=:regular)
        end
            if do_label
                for i in cell_anno
                    anno_df3 = filter(anno => ==(i), anno_df)
                    x_ax = anno_df3[!, x_col]
                    y_ax = anno_df3[!, y_col]
                    MK.text!(ax1, i, position = (mean(x_ax) - label_offset[1], mean(y_ax) - label_offset[2]),align = (:center, :center),font = "Noto Sans Regular",fontsize = label_size,color = label_color)
                end
            end
        else
            MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in polygons]; strokecolor=stroke_color, color=plt_color1, strokewidth=stroke_width)
            if do_legend
                cells = string.(collect(keys(cell_color)))
                colors = collect(values(cell_color))
                for (cell1, color1) in zip(cells, colors)
                    MK.scatter!(ax1, [NaN], [NaN], color = color1, strokewidth = 0.5,strokecolor=stroke_color, markersize = marker_size, label = cell1)
                end
                c_map2 = [(i, alpha) for i in c_map]
                MK.Legend(fig[1, 2], ax1, String(anno), framecolor=:white, labelsize=legend_fontsize, 
                    nbanks=legend_ncol, titlesize=20, titlefont=:regular)
                MK.colgap!(fig.layout, 1)
            end
            MK.colsize!(fig.layout, 1, MK.Aspect(1, 1))
        end

        MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
        MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
        return MK.current_figure()

    elseif data_use == "individual"
        ### xenium processing
        xn_obj = sp.pairedData.xnObj
        if xn_cell_shape == "point"
            img2, anno_df = process_xn_dimplot_data(sp; anno=xn_anno, anno_color=xn_anno_color, x_col = x_col,  y_col = y_col, 
                cell_highlight=xn_cell_highlight, x_lims = x_lims, y_lims = y_lims, pt_bg_color = pt_bg_color, alpha=alpha,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = xn_cell_shape
            )
        else
            img2, polygons, cell_color, plt_color1, c_map = process_xn_dimplot_data(sp; anno=xn_anno, anno_color=xn_anno_color, x_col = x_col,  y_col = y_col, 
                cell_highlight=xn_cell_highlight, x_lims = x_lims, y_lims = y_lims, pt_bg_color = pt_bg_color, alpha=alpha,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = xn_cell_shape
            )
        end
        # visium processing
        if hd_layer == "2_um"
            error("""Your bin size in hd_layer was set to "2_um". Please set it back to "8_um" or "16_um".""")
        end
        sp.pairedData.vsObj = set_default_layer(sp.pairedData.vsObj; layer_slot = hd_layer)
        hd_obj = sp.pairedData.vsObj
        img_vs, poly, cell_color, plt_color = process_hd_dimplot_data(hd_obj; anno=vs_anno, anno_color=vs_anno_color, x_col = x_col, y_col = y_col, pt_bg_color=pt_bg_color, 
            cell_highlight=vs_cell_highlight, x_lims = x_lims, y_lims = y_lims,alpha = alpha, adjust_contrast = adjust_contrast, adjust_brightness = adjust_brightness)
        plt_color=[(i, alpha) for i in plt_color]
        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false);
        ax2 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false)
        MK.Label(fig[0, 1], "Xenium " * xn_anno , fontsize=20, halign=:center, valign=:bottom)
        if plot_img
            MK.image!(ax1, img2)
        end
        if xn_cell_shape == "point"
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
                MK.Legend(fig[1, 2], ax2, String(xn_anno), framecolor=:white, labelsize=legend_fontsize, nbanks=xn_legend_ncol, titlesize=20, titlefont=:regular)
                ax3 = MK.Axis(fig[1,3]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                    xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
                MK.Label(fig[0, 3], "VisiumHD " * vs_anno , fontsize=20, halign=:center, valign=:bottom)
                MK.colsize!(fig.layout, 3, MK.Aspect(1, 1))
            else
                ax3 = MK.Axis(fig[1,2]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                    xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
                MK.Label(fig[0, 2], "VisiumHD " * vs_anno , fontsize=20, halign=:center, valign=:bottom)
                MK.colsize!(fig.layout, 2, MK.Aspect(1, 1))
            end
            if do_label
                for i in cell_anno
                    anno_df3 = filter(anno => ==(i), anno_df)
                    x_ax = anno_df3[!, x_col]
                    y_ax = anno_df3[!, y_col]
                    MK.text!(ax1, i, position = (mean(x_ax) - label_offset[1], mean(y_ax) - label_offset[2]),align = (:center, :center),font = "Noto Sans Regular",fontsize = label_size,color = label_color)
                end
            end
        else
            MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in polygons]; strokecolor=stroke_color, color=plt_color1, strokewidth=stroke_width)
            if do_legend
                cells = string.(collect(keys(cell_color)))
                colors = collect(values(cell_color))
                for (cell1, color1) in zip(cells, colors)
                    MK.scatter!(ax1, [NaN], [NaN], color = color1, strokewidth = 0.5,strokecolor=stroke_color, markersize = marker_size, label = cell1)
                end
                c_map2 = [(i, alpha) for i in c_map]
                MK.Legend(fig[1, 2], ax1, String.(xn_anno), framecolor=:white, labelsize=legend_fontsize, 
                    nbanks=xn_legend_ncol, titlesize=20, titlefont=:regular)
                MK.colgap!(fig.layout, 1)
            end
            MK.colsize!(fig.layout, 1, MK.Aspect(1, 1))
        end
        if do_legend
            ax3 = MK.Axis(fig[1,3]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
            MK.Label(fig[0, 3], "VisiumHD " * vs_anno , fontsize=20, halign=:center, valign=:bottom)
        else
            ax3 = MK.Axis(fig[1,2]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
            MK.Label(fig[0, 2], "VisiumHD " * vs_anno , fontsize=20, halign=:center, valign=:bottom)
        end
        if plot_img
            MK.image!(ax3, img_vs)
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
            MK.Legend(fig[1, 4], ax3, String.(vs_anno), framecolor=:white, labelsize=legend_fontsize, nbanks=vs_legend_ncol, titlesize=20, titlefont=:regular)
            MK.rowgap!(fig.layout, 3)
            MK.colgap!(fig.layout, 1)
            MK.colsize!(fig.layout, 1, MK.Aspect(1, 1))
            MK.colsize!(fig.layout, 3, MK.Aspect(1, 1))
        else
            MK.rowgap!(fig.layout, 3)
            MK.colsize!(fig.layout, 1, MK.Aspect(1, 1))
            MK.colsize!(fig.layout, 2, MK.Aspect(1, 1))
        end
        MK.poly!(ax3, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
        MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
        MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
        return MK.current_figure()
    else
        error("data_use can only be cellseg or individual.")
    end
end

function paired_dim_plot(sp::PairedObject;
    mode::String = "cell",
    bg_mode::String = "gene",
    anno::Union{Symbol, String}="cluster", 
    vs_anno::Union{Symbol, String}="cluster",
    bg_gene::Union{String, Vector{String}, Tuple{String}, Nothing}=nothing,
    anno_color::Union{Nothing, Vector{String}, Tuple{String}} = nothing, 
    vs_anno_color::Union{Nothing, Vector{String}, Tuple{String}} = nothing, 
    color_keys::Union{Vector{String}, Tuple{String}}=["gray94","lemonchiffon","orange","red3"],
    cell_order::Union{Vector{String}, Nothing}=nothing, 
    vs_cell_order::Union{Vector{String}, Nothing}=nothing, 
    cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing, 
    vs_cell_highlight::Union{String, Int64, Vector, Tuple, Nothing}=nothing,
    pt_bg_color="transparent",
    stroke_width=0.5,
    stroke_color=:transparent, 
    label_size=50, 
    label_color="black", 
    do_label=false, 
    alpha::Real = 0.5, 
    legend_ncol = 1,
    plot_img = true,
    img_use::String = "vs_img",
    clip = 0.0,
    scale =  false,
    x_col = "x", 
    y_col = "y",
    hd_layer = "8_um",
    x_lims = nothing, 
    y_lims = nothing,
    marker_size = 2, 
    bg_color = :white,
    gene_colors = nothing,
    bg_tx = false,
    adjust_contrast = 1.0,
    adjust_brightness = 0.0, 
    legend_size = 10, 
    legend_fontsize = 20, 
    do_legend = false,
    height = 500, 
    width = 500
)
    if mode == "cell"
        ## Cell seg data processing
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
            cell_color=Dict(cell_highlight .=> anno_color)
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

        ### background data processing
        hd_obj = sp.pairedData.vsObj
        img, poly2, gene_expr, plt_color, c_map = process_hd_featureplot_data(hd_obj, bg_gene; color_keys = color_keys, x_col = x_col,  
                y_col = y_col, hd_layer = hd_layer, clip = clip,  x_lims = x_lims,  y_lims = y_lims,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness)
        if img_use == "xn_img"
            img = deepcopy(sp.pairedData.xnObj.imageData)
            if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
                img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
            end
            img = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        end
        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, 
            xgridvisible = false, ygridvisible = false)
        if plot_img
            MK.image!(ax1, img)
        end
        MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in polygons]; strokecolor=stroke_color, color=plt_color1, strokewidth=stroke_width)
        MK.Label(fig[0, 1], "VisiumHD " * bg_gene * " transcript on Xenium cell seg", fontsize=18, halign=:center, valign=:bottom)
        if do_legend
            cells = string.(collect(keys(cell_color)))
            colors = collect(values(cell_color))
            for (cell1, color1) in zip(cells, colors)
                MK.scatter!(ax1, [NaN], [NaN], color = color1, strokewidth = 0.5,strokecolor=stroke_color, markersize = 3*legend_size, label = cell1)
            end
            c_map2 = [(i, alpha) for i in c_map]
            MK.Colorbar(fig[1,2], colormap = c_map2,  width=15, limits = (0, maximum(gene_expr)))
            MK.Label(fig[0, 2], bg_gene, fontsize=18)
            MK.Legend(fig[1, 3], ax1, String.(anno), framecolor=:white, labelsize=legend_fontsize, 
                nbanks=legend_ncol, titlesize=20, titlefont=:regular)
            MK.colgap!(fig.layout, 1)
        end
        MK.colsize!(fig.layout, 1, MK.Aspect(1, 1.1))
        MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in poly2]; strokecolor=stroke_color, 
                color=plt_color, strokewidth=stroke_width,label="")
        MK.rowgap!(fig.layout, 3)
        MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
        MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
        return MK.current_figure() 
    elseif mode == "bin"
        if hd_layer == "2_um"
            error("""Your bin size in hd_layer was set to "2_um". Please set it back to "8_um" or "16_um".""")
        end
        sp.pairedData.vsObj = set_default_layer(sp.pairedData.vsObj; layer_slot = hd_layer)
        hd_obj = sp.pairedData.vsObj
        img, poly, cell_color, plt_color = process_hd_dimplot_data(hd_obj; anno=vs_anno, anno_color=vs_anno_color, x_col = x_col, y_col = y_col, pt_bg_color=pt_bg_color, 
            cell_highlight=vs_cell_highlight, x_lims = x_lims, y_lims = y_lims,alpha = alpha, adjust_contrast = adjust_contrast, adjust_brightness = adjust_brightness)
        if img_use == "xn_img"
            img = deepcopy(sp.pairedData.xnObj.imageData)
            if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
                img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
            end
            img = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        end
        xn_obj = sp.pairedData.xnObj
        df_plt, all_genes, all_colors = process_xn_transcript_data(xn_obj, bg_gene; x_lims = x_lims, y_lims = y_lims, x_col = x_col,  
                            y_col = y_col,  gene_colors = gene_colors, bg_tx = bg_tx)
        all_colors=[(i, alpha) for i in all_colors]
        #plt_color=[(i, alpha) for i in plt_color]
        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false)
        ax2 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false)
        if plot_img
            MK.image!(ax1, img)
        end
        MK.Label(fig[0, 1], "Xenium transcript(s) on VisiumHD cell bin", fontsize=18, halign=:center, valign=:bottom)
        if do_legend
            if !isa(vs_cell_order, Nothing)
            cells = vs_cell_order
            else
            cells = string.(collect(keys(cell_color)))
            end
            colors = [cell_color[i] for i in cells]
            for (cell1, color1) in zip(cells, colors)
                MK.scatter!(ax1,[NaN], [NaN], color = color1, marker=:rect,
                                strokewidth = 0.5,strokecolor=stroke_color, markersize = 2*marker_size, label = cell1)
            end
            MK.Legend(fig[1, 2], ax1, vs_anno, framecolor=:white, labelsize=legend_fontsize, nbanks=legend_ncol, titlesize=20, titlefont=:regular)
        end
        MK.colsize!(fig.layout, 1, MK.Aspect(1, 1.1))
        MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
        for (gene, ann_color) in zip(all_genes, all_colors)
            x_ax = df_plt[!, x_col][df_plt.new_gene .== gene]
            y_ax = df_plt[!, y_col][df_plt.new_gene .== gene]
            if do_legend
                MK.scatter!(ax2, x_ax , y_ax;  visible=false,
                                color=ann_color, strokewidth=0, markersize=3*legend_size, label=gene)
                MK.scatter!(ax1, x_ax , y_ax; color = ann_color, strokewidth = 0, markersize = marker_size)
            else
                MK.scatter!(ax1, x_ax , y_ax; color = ann_color, strokewidth = 0, markersize = marker_size)
            end
        end
        if do_legend
            MK.Legend(fig[1, 3], ax2, "Transcript", framecolor=:white, labelsize=legend_fontsize, titlesize=20, titlefont=:regular)
        end
        MK.rowgap!(fig.layout, 3)
        MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
        MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
        return MK.current_figure()
    else
        error("""The parameter mode can only be "cell" or "bin"!""")
    end
end

function paired_feature_plot(sp::PairedObject, gene::String;
    mode = "bin",
    tx_list::Union{String, Vector{String}, Tuple{String}, Nothing}=nothing,
    color_keys::Union{Vector{String}, Tuple{String}}=["gray94","lemonchiffon","orange","red3"],
    x_col = "x",  
    y_col = "y", 
    hd_layer = "8_um",
    clip = 0.0,
    scale =  false,
    x_lims = nothing, 
    y_lims = nothing,
    adjust_contrast= 1.0,
    adjust_brightness = 0.0,
    img_use = "xn_img",
    plot_img = true,
    marker_size = 10, 
    bg_color = :white,
    gene_colors = nothing,
    bg_tx = false,
    legend_size = 10, 
    legend_fontsize = 20, 
    do_legend = false,
    alpha::Real = 0.5, 
    stroke_width=0.5,
    stroke_color=:transparent, 
    height = 500, 
    width = 500        
)
    if mode == "bin"
        hd_obj = sp.pairedData.vsObj
        img, poly, gene_expr, plt_color, c_map = process_hd_featureplot_data(hd_obj, gene; color_keys = color_keys, x_col = x_col,  
                y_col = y_col, hd_layer = hd_layer, clip = clip,  x_lims = x_lims,  y_lims = y_lims,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness)
        if img_use == "xn_img"
            img = deepcopy(sp.pairedData.xnObj.imageData)
            if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
                img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
            end
            img = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        end
    elseif mode == "cell"
        img, poly, gene_expr, plt_color, c_map = process_paired_featureplot_data(sp, gene; color_keys = color_keys, x_col = x_col,  
            y_col = y_col, clip = clip,  x_lims = x_lims,  y_lims = y_lims, cell_shape = "point",
            adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, img_use = img_use)
    else
        error("""The parameter mode can only be "cell" or "bin"!""")
    end
    if isa(tx_list, Nothing)
        tx_list = gene
    end
    xn_obj = sp.pairedData.xnObj
    df_plt, all_genes, all_colors = process_xn_transcript_data(xn_obj, tx_list; x_lims = x_lims, y_lims = y_lims, x_col = x_col,  
                    y_col = y_col,  gene_colors = gene_colors, bg_tx = bg_tx)
    all_colors=[(i, alpha) for i in all_colors]
    plt_color=[(i, alpha) for i in plt_color]
    fig = MK.Figure(size=(width, height))
    ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, 
        xgridvisible = false, ygridvisible = false)
    ax2 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
        xgridvisible = false,ygridvisible = false)        
    if plot_img
        MK.image!(ax1, img)
    end
    MK.Label(fig[0, 1], "Xenium transcript(s) on " * "VisiumHD " * gene * " expression" , fontsize=18, halign=:center, valign=:bottom)
    MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, 
            color=plt_color, strokewidth=stroke_width,label="")
    if do_legend
        c_map2 = [(i, alpha) for i in c_map]
        MK.Colorbar(fig[1,2], colormap = c_map2,  width=15, limits = (0, maximum(gene_expr)))
        MK.Label(fig[0, 2], gene, fontsize=18)
        MK.colgap!(fig.layout, 1)
    end
    MK.colsize!(fig.layout, 1, MK.Aspect(1, 1.1))
    for (gene, ann_color) in zip(all_genes, all_colors)
        x_ax = df_plt[!, x_col][df_plt.new_gene .== gene]
        y_ax = df_plt[!, y_col][df_plt.new_gene .== gene]
        if do_legend
            MK.scatter!(ax2, x_ax , y_ax;  visible=false,
                            color=ann_color, strokewidth=0, markersize=2*legend_size, label=gene)
            MK.scatter!(ax1, x_ax , y_ax; color = ann_color, strokewidth = 0, markersize = marker_size)
        else
            MK.scatter!(ax1, x_ax , y_ax; color = ann_color, strokewidth = 0, markersize = marker_size)
        end
    end
    if do_legend
        MK.Legend(fig[1, 3], ax2, "Transcript", framecolor=:white, labelsize=legend_fontsize, titlesize=20, titlefont=:regular)
    end
    MK.rowgap!(fig.layout, 3)
    MK.xlims!(MK.current_axis(), x_lims .- x_lims[1] .+ 1)
    MK.ylims!(MK.current_axis(), y_lims .- y_lims[1] .+ 1)
    return MK.current_figure()            
end

function gemini_dim_plot(sp::PairedObject; 
    xn_anno::Union{Symbol, String}="cluster", 
    vs_anno::Union{Symbol, String}="cluster", 
    xn_anno_color::Union{Vector{String}, Nothing}= nothing, 
    vs_anno_color::Union{Vector{String}, Nothing} = nothing, 
    xn_cell_order::Union{Vector{String}, Nothing}=nothing, 
    vs_cell_order::Union{Vector{String}, Nothing}=nothing,
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
    hd_layer = "8_um",
    marker_size = 2, 
    bg_color = :white,  
    adjust_contrast = 1.0,
    adjust_brightness = 0.0, 
    legend_size = 30, 
    legend_fontsize = 20, 
    do_legend = false,
    height = 500, 
    width = 500,
    break_ratio = 0.5,
    aspect_ratio = 0.8,
    xn_cell_shape = "point"
)
        if break_ratio < 0.05 || break_ratio > 0.95
            error("break_ratio should not be < 0.05 or > 0.95")
        end
        x_coord_xn = deepcopy(sp.spmetaData.cell.x)
        y_coord_xn = deepcopy(sp.spmetaData.cell.y)
        x_coord_vs = deepcopy(sp.pairedData.vsObj.spmetaData.pxl_row_in_fullres)
        y_coord_vs = deepcopy(sp.pairedData.vsObj.spmetaData.pxl_col_in_fullres)
        x_lims_xn=(minimum(x_coord_xn)-0.05*maximum(x_coord_xn), maximum(x_coord_xn) * break_ratio)
        x_lims_xn=adjust_lims(x_lims_xn)
        x_lims_vs=(maximum(x_coord_vs) * break_ratio, size(sp.pairedData.vsObj.imageData.fullresImage)[1])
        x_lims_vs=adjust_lims(x_lims_vs)
        y_lims=(minimum(y_coord_xn)-0.05*maximum(y_coord_xn),1.05*maximum(y_coord_xn))
        y_lims2=(minimum(y_coord_vs)-0.05*maximum(y_coord_vs),1.05*maximum(y_coord_vs))
        if xn_cell_shape == "point"
            img2, anno_df = process_xn_dimplot_data(sp; anno=xn_anno, anno_color=xn_anno_color, x_col = x_col,  y_col = y_col, 
                cell_highlight=xn_cell_highlight, x_lims = x_lims_xn, y_lims = y_lims, pt_bg_color = pt_bg_color, alpha=alpha,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = xn_cell_shape
            )
        else
            img2, polygons, cell_color, plt_color1, c_map = process_xn_dimplot_data(sp; anno=xn_anno, anno_color=xn_anno_color, x_col = x_col,  y_col = y_col, 
                cell_highlight=xn_cell_highlight, x_lims = x_lims_xn, y_lims = y_lims, pt_bg_color = pt_bg_color, alpha=alpha,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = xn_cell_shape
            )
        end
        img2 = flip_bg_color(img2)
        if hd_layer == "2_um"
            error("""Your bin size in hd_layer was set to "2_um". Please set it back to "8_um" or "16_um".""")
        end
        sp.pairedData.vsObj = set_default_layer(sp.pairedData.vsObj; layer_slot = hd_layer)
        hd_obj = sp.pairedData.vsObj
        img_vs, poly, cell_color, plt_color = process_hd_dimplot_data(hd_obj; anno=vs_anno, anno_color=vs_anno_color, x_col = x_col, y_col = y_col, pt_bg_color=pt_bg_color, 
            cell_highlight=vs_cell_highlight, x_lims = x_lims_vs, y_lims = y_lims,alpha = alpha, adjust_contrast = adjust_contrast, adjust_brightness = adjust_brightness)
        plt_color=[(i, alpha) for i in plt_color]
        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, 
            xautolimitmargin = (0.05, 0),
            xgridvisible = false,ygridvisible = false)
        ax2 = MK.Axis(fig[1,1]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xautolimitmargin = (0.05, 0), 
            xgridvisible = false,ygridvisible = false)
        MK.Label(fig[0, 1], "Xenium", fontsize=20, halign=:center, valign=:bottom)
        if plot_img
            MK.image!(ax1, img2)
        end
        if xn_cell_shape == "point"
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
                MK.Legend(fig[1, 0],ax2, String(xn_anno), framecolor=:white, labelsize=legend_fontsize, nbanks=xn_legend_ncol, titlesize=20, titlefont=:regular)
            end
            if do_label
                for i in cell_anno
                    anno_df3 = filter(anno => ==(i), anno_df)
                    x_ax = anno_df3[!, x_col]
                    y_ax = anno_df3[!, y_col]
                    MK.text!(ax1, i, position = (mean(x_ax) - label_offset[1], mean(y_ax) - label_offset[2]),align = (:center, :center),font = "Noto Sans Regular",fontsize = label_size,color = label_color)
                end
            end
        else
            MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in polygons]; strokecolor=stroke_color, color=plt_color1, strokewidth=stroke_width)
            if do_legend
                cells = string.(collect(keys(cell_color)))
                colors = collect(values(cell_color))
                for (cell1, color1) in zip(cells, colors)
                    MK.scatter!(ax1, [NaN], [NaN], color = color1, strokewidth = 0.5,strokecolor=stroke_color, markersize = marker_size, label = cell1)
                end
                c_map2 = [(i, alpha) for i in c_map]
                MK.Legend(fig[1, 0], ax1, String.(xn_anno), framecolor=:white, labelsize=legend_fontsize, 
                    nbanks=xn_legend_ncol, titlesize=20, titlefont=:regular)
                MK.colgap!(fig.layout, 1)
            end
            MK.colsize!(fig.layout, 1, MK.Aspect(1, aspect_ratio))
        end
        ax3 = MK.Axis(fig[1,2]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xautolimitmargin = (0, 0.05), 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
        MK.Label(fig[0, 2], "VisiumHD", fontsize=20, halign=:center, valign=:bottom)
        MK.colsize!(fig.layout, 2, MK.Aspect(1, aspect_ratio))
        if plot_img
            MK.image!(ax3, img_vs)
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
            MK.Legend(fig[1, 3], ax3, String.(vs_anno), framecolor=:white, labelsize=legend_fontsize, nbanks=vs_legend_ncol, titlesize=20, titlefont=:regular)
            MK.rowgap!(fig.layout, 2)
            MK.colgap!(fig.layout, 0)
            MK.colsize!(fig.layout, 1, MK.Aspect(1, aspect_ratio))
            MK.colsize!(fig.layout, 2, MK.Aspect(1, aspect_ratio))
        else
            MK.rowgap!(fig.layout, 2)
            MK.colgap!(fig.layout, 0)
            MK.colsize!(fig.layout, 1, MK.Aspect(1, aspect_ratio))
            MK.colsize!(fig.layout, 2, MK.Aspect(1, aspect_ratio))
        end
        MK.poly!(ax3, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
        y_lims = [max(y_lims[1], y_lims2[1]), min(y_lims[2], y_lims2[2])]
        y_lims[1] = y_lims[1] > 0 ? 0 : y_lims[1]
        MK.ylims!(ax1, y_lims...)
        MK.ylims!(ax3, y_lims...)
        return MK.current_figure()
end