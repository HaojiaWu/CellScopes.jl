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
            cell_highlight=vs_cell_highlight, x_lims = x_lims, y_lims = y_lims,alpha = alpha, adjust_contrast = adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = "bin")
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
        return fig
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
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = "bin")
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
            cell_highlight=vs_cell_highlight, x_lims = x_lims, y_lims = y_lims,alpha = alpha, adjust_contrast = adjust_contrast, adjust_brightness = adjust_brightness, cell_shape="bin")
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
        return fig
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
    x_lims = nothing, 
    y_lims = nothing,
    adjust_contrast= 1.0,
    adjust_brightness = 0.0,
    img_use = "xn_img",
    cell_shape ="point",
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

    if gene ∉ sp.rawCount.gene_name
        error("Gene not found: $gene. This gene may not be detected in this dataset.")
    end
    if mode == "bin"
        hd_obj = sp.pairedData.vsObj
        img, poly, gene_expr, plt_color, c_map = process_hd_featureplot_data(hd_obj, gene; color_keys = color_keys, x_col = x_col,  
                y_col = y_col, hd_layer = hd_layer, clip = clip,  x_lims = x_lims,  y_lims = y_lims,
                adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = "bin")
        if img_use == "xn_img"
            img = deepcopy(sp.pairedData.xnObj.imageData)
            if !isa(x_lims, Nothing) && !isa(y_lims, Nothing)
                img = img[round(Int,x_lims[1]):round(Int, x_lims[2]), round(Int, y_lims[1]):round(Int, y_lims[2])]
            end
            img = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        end
    elseif mode == "cell"
        img, poly, gene_expr, plt_color, c_map = process_paired_featureplot_data(sp, gene; color_keys = color_keys, x_col = x_col,  
            y_col = y_col, clip = clip,  x_lims = x_lims,  y_lims = y_lims, cell_shape = cell_shape,
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
    return fig            
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
    margin = 0.01,
    hd_layer = "8_um",
    only_selected = true,
    marker_size = 2, 
    canvas_color= :white,
    bg_color = :gray85,  
    adjust_contrast = 1.0,
    adjust_brightness = 0.0, 
    legend_size = 30, 
    legend_fontsize = 20, 
    do_legend = false,
    height = 1000, 
    width = 600,
    break_ratio = 0.5,
    aspect_ratio = 0.8,
    cell_shape = "point"
)
        if break_ratio < 0.05 || break_ratio > 0.95
            error("break_ratio should not be < 0.05 or > 0.95")
        end
        x_coord_xn = deepcopy(sp.spmetaData.cell.x)
        y_coord_xn = deepcopy(sp.spmetaData.cell.y)
        x_coord_vs = deepcopy(sp.pairedData.vsObj.spmetaData.pxl_row_in_fullres)
        y_coord_vs = deepcopy(sp.pairedData.vsObj.spmetaData.pxl_col_in_fullres)
        img_vs_size = size(sp.pairedData.vsObj.imageData.fullresImage)
        img_xn_size = size(sp.pairedData.xnObj.imageData)
        img_limit = [maximum([img_vs_size[1], img_xn_size[1]]), minimum([img_vs_size[2], img_xn_size[2]])]
        x_lims_xn = [minimum(x_coord_xn)-margin*maximum(x_coord_xn), maximum(x_coord_xn) * break_ratio]
        if x_lims_xn[1] < 1
            x_lims_xn[1] = 1
        end
        x_lims_xn=adjust_lims(x_lims_xn)
        x_lims_vs=[maximum(x_coord_xn) * break_ratio, img_limit[1]-1]
        x_lims_vs=adjust_lims(x_lims_vs)
        y_lims=[minimum(y_coord_xn)-margin*maximum(y_coord_xn),(1.0+margin)*maximum(y_coord_xn)]
        if y_lims[1] < 1
            y_lims[1] = 1
        end
        y_lims[2] = y_lims[2] > img_limit[2] ? (img_limit[2] -1) : y_lims[2]
        y_lims2=[minimum(y_coord_vs)-margin*maximum(y_coord_vs),(1.0+margin)*maximum(y_coord_vs)]
        if y_lims2[1] < 1
            y_lims2[1] = 1
        end
        y_lims2[2] = y_lims2[2] > img_limit[2] ? (img_limit[2] -1) : y_lims2[2]
 
        img2, anno_df = process_xn_dimplot_data(sp; anno=xn_anno, anno_color=xn_anno_color, x_col = x_col,  y_col = y_col, 
            cell_highlight=xn_cell_highlight, x_lims = x_lims_xn, y_lims = y_lims, pt_bg_color = pt_bg_color, alpha=alpha,
            adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, cell_shape = "point"
        )

        flip_bg_color!(img2)
        if hd_layer == "2_um"
            error("""Your bin size in hd_layer was set to "2_um". Please set it back to "8_um" or "16_um".""")
        end
        sp.pairedData.vsObj = set_default_layer(sp.pairedData.vsObj; layer_slot = hd_layer)
        hd_obj = sp.pairedData.vsObj
        img_vs, poly, cell_color, plt_color = process_hd_dimplot_data(hd_obj; anno=vs_anno, anno_color=vs_anno_color, x_col = x_col, y_col = y_col, pt_bg_color=pt_bg_color, 
            cell_highlight=vs_cell_highlight, x_lims = x_lims_vs, y_lims = y_lims,alpha = alpha, adjust_contrast = adjust_contrast, adjust_brightness = adjust_brightness,
            cell_shape = cell_shape)
        plt_color=[(i, alpha) for i in plt_color]
        fig = MK.Figure(size=(width, height))
        ax1 = MK.Axis(fig[1,1]; backgroundcolor = canvas_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, 
            xautolimitmargin = (0.05, 0),
            xgridvisible = false,ygridvisible = false)
        ax2 = MK.Axis(fig[1,1]; backgroundcolor = canvas_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xautolimitmargin = (0.05, 0), 
            xgridvisible = false,ygridvisible = false)
        MK.Label(fig[0, 1], "Xenium", fontsize=20, halign=:center, valign=:bottom)
        if plot_img
            MK.image!(ax1, img2)
        end
        if isa(xn_cell_order, Nothing)
            cell_anno=unique(anno_df[!,xn_anno])
        else
            cell_anno = xn_cell_order
        end
        if !only_selected
            bg_cells = deepcopy(sp.spmetaData.cell)
            bg_cells = filter(:cell => !(∈(Set(anno_df.cell))), bg_cells)
            bg_cells = filter([:x, :y] => (x, y) -> x_lims_xn[1] < x < x_lims_xn[2] && y_lims[1] < y < y_lims[2], bg_cells)
            bg_cells[!, x_col] = bg_cells[!, x_col] .- x_lims_xn[1]
            bg_cells[!, y_col] = bg_cells[!, y_col] .- y_lims[1]
            MK.scatter!(ax1, bg_cells[!, x_col], bg_cells[!, y_col]; color = bg_color, strokewidth = 0, markersize = marker_size)
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

        ax3 = MK.Axis(fig[1,2]; backgroundcolor = canvas_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
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
    if !only_selected
        anno_df = deepcopy(hd_obj.spmetaData)
        all_poly = deepcopy(hd_obj.polygonData)
        all_poly = [m .- [x_lims_vs[1]-1 y_lims[1]-1] for m in all_poly]
        rename!(anno_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, :x, :y])
        select_fov = deepcopy(anno_df)
        select_fov = filter(vs_anno => !(∈(Set(vs_cell_highlight))), select_fov)
        select_fov = filter([:x, :y] => (x, y) -> x_lims_vs[1] < x < x_lims_vs[2] && y_lims[1] < y < y_lims[2], select_fov)
        select_fov[!, x_col] = select_fov[!, x_col] .- x_lims_vs[1]
        select_fov[!, y_col] = select_fov[!, y_col] .- y_lims[1]
        if cell_shape != "point"
            polygon_num = select_fov.ID
            bg_poly = all_poly[polygon_num]
            MK.poly!(ax3, [MK.Point2.(eachrow(p)) for p in bg_poly]; strokecolor=stroke_color, 
                    color=bg_color, strokewidth=stroke_width,label="")
        else
            MK.scatter!(ax3, select_fov[!, x_col], select_fov[!, y_col]; color = bg_color, strokewidth = 0, markersize = marker_size)
        end
    end
    if cell_shape !="point"
        MK.poly!(ax3, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
    else
        poly[!, x_col] = poly[!, x_col] .- x_lims_vs[1]
        poly[!, y_col] = poly[!, y_col] .- y_lims[1]
        MK.scatter!(ax3, poly[!, x_col], poly[!, y_col]; color = poly.new_color, strokewidth = 0, markersize = marker_size)
    end
    y_lims = [max(y_lims[1], y_lims2[1]), min(y_lims[2], y_lims2[2])]
    y_lims[1] = y_lims[1] > 0 ? 0 : y_lims[1]
    MK.ylims!(ax1, y_lims...)
    MK.ylims!(ax3, y_lims...)
    return fig
end

function gemini_feature_plot(sp, gene::String;
    color_keys_xn::Union{Vector{String}, Tuple{String}}=["gray85","cyan","blue","blue3"],
    color_keys_vs::Union{Vector{String}, Tuple{String}}=["gray85","green3","darkgreen","#013300"],
    x_col = "x",  
    y_col = "y", 
    hd_layer = "8_um",
    clip = 0.0,
    order=true,
    adjust_contrast= 1.0,
    adjust_brightness = 0.0,
    plot_img = true,
    marker_size = 2, 
    canvas_color = :white,
    bg_color = :gray85,
    do_legend = false,
    alpha::Real = 0.5, 
    stroke_width=0.5,
    cell_shape = "point",
    margin = 0.01,
    stroke_color=:transparent,
    break_ratio=0.5,
    aspect_ratio=0.7,
    only_expr=false,
    horizontal_cut = false,
    height = 1000, 
    width = 600        
)
    if gene ∉ get_shared_gene(sp)
        error("The gene '$gene' is not found among the shared genes. It may not be commonly detected across the datasets in this object.")
    end
    if break_ratio < 0.05 || break_ratio > 0.95
        error("break_ratio should not be < 0.05 or > 0.95")
    end
    img_xn_size = size(sp.pairedData.xnObj.imageData)
    img_vs_size = size(sp.pairedData.vsObj.imageData.fullresImage)
    img_limit = [minimum([img_vs_size[1], img_xn_size[1]]), minimum([img_vs_size[2], img_xn_size[2]])]
    if !horizontal_cut
        x_lims_xn = [1, img_limit[1] * break_ratio]
        x_lims_xn = adjust_lims(x_lims_xn)
        x_lims_vs = [img_limit[1] * break_ratio, img_limit[1]-1]
        x_lims_vs = adjust_lims(x_lims_vs)
        y_lims_vs = y_lims_xn = [1, img_limit[2] - 1]
    else
        y_lims_vs = [1, img_limit[2] * break_ratio]
        y_lims_vs = adjust_lims(y_lims_vs)
        y_lims_xn = [img_limit[2] * break_ratio, img_limit[2]-1]
        y_lims_xn = adjust_lims(y_lims_xn)
        x_lims_vs = x_lims_xn = [1, img_limit[1] - 1]
    end

    # visiumHD gene expr processing
    hd_obj = sp.pairedData.vsObj
    img, poly, gene_expr, plt_color, c_map = process_hd_featureplot_data(hd_obj, gene; color_keys = color_keys_vs, x_col = x_col,  
            y_col = y_col, hd_layer = hd_layer, clip = clip,  x_lims = x_lims_vs,  y_lims = y_lims_vs, cell_shape = cell_shape,
            adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness)

    # xenium gene expr processing
    xn_obj = sp.pairedData.xnObj
    img_xn, df_plt, gene_expr_xn, plt_color_xn, c_map_xn = process_paired_featureplot_data(xn_obj, gene; color_keys = color_keys_xn, x_col = x_col,  
        y_col = y_col, clip = clip,  x_lims = x_lims_xn,  y_lims = y_lims_xn, cell_shape = "point",
        adjust_contrast= adjust_contrast, adjust_brightness = adjust_brightness, img_use = "xn_img")
    flip_bg_color!(img_xn)
    plt_color_xn=[(i, alpha) for i in plt_color_xn]
    df_plt.gene .= gene_expr_xn
    if sum(gene_expr) > 0.0
        df_plt.plt_color = plt_color_xn
        if order
            
            df_plt = sort(df_plt,:gene);
        end
    else
        plt_color_xn = repeat([color_keys_xn[1]], length(gene_expr_xn))
        df_plt.plt_color = plt_color_xn
    end
    df_plt[!, x_col] = df_plt[!, x_col] .- x_lims_xn[1]
    df_plt[!, y_col] = df_plt[!, y_col] .- y_lims_xn[1]

    fig = MK.Figure(size=(width, height))
    if !horizontal_cut
        ax1_pos = fig[1, 1]; label1_pos = fig[0, 1]; cbar1_pos = fig[1, 0]; gene_label1_pos = fig[0, 0]
        ax2_pos = fig[1, 2]; label2_pos = fig[0, 2]; cbar2_pos = fig[1, 3]; gene_label2_pos = fig[0, 3]
        label_rot = 0; cbar_horizontal = true; label_halign = :center; label_valign = :bottom; cbar_width = 10
        cbar_tickalign_xn = 0; cbar_flipaxis_xn = false; cbar_tickalign_vs = 1; cbar_flipaxis_vs = true
    else
        ax1_pos = fig[1, 1]; label1_pos = fig[1, 0]; cbar1_pos = fig[0, 1]; gene_label1_pos = fig[0, 0]
        ax2_pos = fig[2, 1]; label2_pos = fig[2, 0]; cbar2_pos = fig[3, 1]; gene_label2_pos = fig[3, 0]
        label_rot = π/2; cbar_horizontal = false; label_halign = :right; label_valign = :center; cbar_width = 0.9 * width
        cbar_tickalign_xn = 0; cbar_flipaxis_xn = true; cbar_tickalign_vs = 0; cbar_flipaxis_vs = false
    end
    ax1 = MK.Axis(ax1_pos; backgroundcolor = canvas_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xautolimitmargin = (0.05, 0), 
        xgridvisible = false, ygridvisible = false)      
    if plot_img
        MK.image!(ax1, img_xn)
    end
    MK.Label(label1_pos, "Xenium", fontsize=18, halign=label_halign, valign=label_valign, rotation=label_rot)
    if !only_expr
        bg_cells = deepcopy(sp.spmetaData.cell)
        bg_cells = filter(:cell => !(∈(Set(df_plt.cell))), bg_cells)
        bg_cells= filter([:x, :y] => (x, y) -> x_lims_xn[1] < x < x_lims_xn[2] && y_lims_xn[1] < y < y_lims_xn[2], bg_cells)
        bg_cells[!, x_col] = bg_cells[!, x_col] .- x_lims_xn[1]
        bg_cells[!, y_col] = bg_cells[!, y_col] .- y_lims_xn[1]
        MK.scatter!(ax1, bg_cells[!, x_col], bg_cells[!, y_col]; color = bg_color, strokewidth = 0, markersize = marker_size)
    end
    MK.scatter!(ax1, df_plt[!, x_col], df_plt[!, y_col]; color = df_plt.plt_color, strokewidth = 0, markersize = marker_size)
    if do_legend
        c_map2 = [(i, alpha) for i in c_map_xn]
        MK.Colorbar(cbar1_pos, colormap = c_map2,  width=cbar_width, limits = (0, maximum(gene_expr_xn)), 
            tickalign = cbar_tickalign_xn, flipaxis = cbar_flipaxis_xn, vertical= cbar_horizontal)
        MK.Label(gene_label1_pos, gene, fontsize=16, halign=label_halign, valign=label_valign, rotation=label_rot)
    end

    ax2 = MK.Axis(ax2_pos; backgroundcolor = canvas_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xautolimitmargin = (0, 0.05), 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false)
    if plot_img
        MK.image!(ax2, img)
    end
    MK.Label(label2_pos, "VisiumHD", fontsize=18, halign=label_halign, valign=label_valign, rotation=label_rot)
    if !only_expr
        anno_df = deepcopy(hd_obj.spmetaData)
        all_poly = deepcopy(hd_obj.polygonData)
        all_poly = [m .- [x_lims_vs[1]-1 y_lims_vs[1]-1] for m in all_poly]
        rename!(anno_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, :x, :y])
        norm_count=hd_obj.normCount
        gene_expr2 = subset_count(norm_count; genes = [gene])
        gene_expr2 = (vec ∘ collect)(gene_expr2.count_mtx)
        anno_df.gene = gene_expr2
        select_fov = anno_df
        select_fov = filter([:x, :y, :gene] => (x, y, gene) -> x_lims_vs[1] < x < x_lims_vs[2] && y_lims_vs[1] < y < y_lims_vs[2] && gene <= clip, select_fov)
        select_fov[!, x_col] = select_fov[!, x_col] .- x_lims_vs[1]
        select_fov[!, y_col] = select_fov[!, y_col] .- y_lims_vs[1]
        if cell_shape != "point"
            polygon_num = select_fov.ID
            bg_poly = all_poly[polygon_num]
            MK.poly!(ax2, [MK.Point2.(eachrow(p)) for p in bg_poly]; strokecolor=stroke_color, 
                    color=bg_color, strokewidth=stroke_width,label="")
        else
            MK.scatter!(ax2, select_fov[!, x_col], select_fov[!, y_col]; color = bg_color, strokewidth = 0, markersize = marker_size)
        end
    end
    plt_color=[(i, alpha) for i in plt_color]
    poly.gene = gene_expr
    if sum(gene_expr) > 0.0
        poly.plt_color = plt_color
        if order
            poly = sort(poly,:gene);
        end
    else
        plt_color = repeat([plt_color[1]], length(gene_expr))
        poly.plt_color = plt_color
    end
    if cell_shape != "point"
        MK.poly!(ax2, [MK.Point2.(eachrow(p)) for p in poly]; strokecolor=stroke_color, 
                color=plt_color, strokewidth=stroke_width,label="")
    else
        poly[!, x_col] = poly[!, x_col] .- x_lims_vs[1]
        poly[!, y_col] = poly[!, y_col] .- y_lims_vs[1]
        MK.scatter!(ax2, poly[!, x_col], poly[!, y_col]; color = poly.plt_color, strokewidth = 0, markersize = marker_size)
    end
    if do_legend
        c_map2 = [(i, alpha) for i in c_map]
        MK.Colorbar(cbar2_pos, colormap = c_map2,  width=cbar_width, limits = (0, maximum(gene_expr)),
            tickalign = cbar_tickalign_vs, flipaxis = cbar_flipaxis_vs, vertical= cbar_horizontal)
        MK.Label(gene_label2_pos, gene, fontsize=16, halign=label_halign, valign=label_valign, rotation=label_rot)
    end
    if !horizontal_cut
        MK.colsize!(fig.layout, 1, MK.Aspect(1, aspect_ratio))
        MK.colsize!(fig.layout, 1, MK.Relative(break_ratio))
        MK.colsize!(fig.layout, 2, MK.Aspect(1, aspect_ratio))
        MK.colsize!(fig.layout, 2, MK.Relative(1 - break_ratio))
        MK.colsize!(fig.layout, 0, MK.Relative(0.02))
        MK.colsize!(fig.layout, 3, MK.Relative(0.02))
    else
        MK.colsize!(fig.layout, 0, MK.Relative(0.05))
        MK.colsize!(fig.layout, 1, MK.Relative(0.95))
        MK.rowsize!(fig.layout, 0, MK.Relative(0.05))
        MK.rowsize!(fig.layout, 1, MK.Relative(0.45))
        MK.rowsize!(fig.layout, 3, MK.Relative(0.05))
        MK.rowsize!(fig.layout, 2, MK.Relative(0.45))

    end
    MK.rowgap!(fig.layout, 0)
    MK.colgap!(fig.layout, 0)
    return fig
end

function sp_feature_plot(sp::IntegratedObject, genes::Union{String, Vector{String}}; count_type = "norm",x_col="x", y_col="y", marker_size=4, order=true,
    color_keys::Union{Vector{String}, Tuple{String,String,String}}=("black","yellow","red"), do_dimname::Bool=false,
        titlesize::Int64 = 24, height::Real = 500, width::Real = 500)
        dim_data = deepcopy(sp.spmetaData)
        if isa(genes, String)
          genes = [genes]
        end
        if count_type === "norm"
            ct_obj = subset_count(sp.normCount; genes = genes)
        elseif count_type === "raw"
            ct_obj = subset_count(sp.rawCount; genes = genes)
        elseif count_type === "scale"
            ct_obj = subset_count(sp.scaleCount; genes = genes)
        else
            println("count_type can only be \"raw\", \"norm\" or \"scale\"!")
        end
        count_mat = ct_obj.count_mtx
        count_mat = DataFrame(count_mat, :auto)
        count_mat.gene = ct_obj.gene_name
        count_mat = permutedims(count_mat, :gene)
        count_mat.cells = sp.rawCount.cell_name
        gene_data = [dim_data count_mat]
        x_lims=(minimum(gene_data[!, x_col])-0.05*maximum(gene_data[!, x_col]),1.05*maximum(gene_data[!, x_col]))
        y_lims=(minimum(gene_data[!, y_col])-0.05*maximum(gene_data[!, y_col]),1.05*maximum(gene_data[!, y_col]))
        c_map = ColorSchemes.ColorScheme([parse(Colorant, color_keys[1]),parse(Colorant, color_keys[2]),parse(Colorant, color_keys[3])])
            group_arr = string.(sp.spmetaData[!, :dataset])
            group_names = unique(group_arr)
            gene_data[!, :dataset] = group_arr
            fig = MK.Figure(size = (width * length(group_names), height * length(genes)))
            for (i, group) in enumerate(group_names)
                for (j, gene) in enumerate(genes)
                    df_plt = gene_data[!, [x_col, y_col, gene, "dataset"]]
                    gene_expr = float.(df_plt[!, gene])
                    if sum(gene_expr) !== 0.0
                        @inbounds colors = get(c_map, gene_expr, :extrema)
                        plt_color = "#" .* hex.(colors)
                        df_plt.plt_color = plt_color
                        if order
                            df_plt = sort(df_plt, Symbol(gene))
                        end
                    else
                        @inbounds plt_color = repeat([color_keys[1]], length(gene_expr))
                        df_plt.plt_color = plt_color
                    end
                    df_plt=filter(:dataset => ==(group), df_plt)
                    if i == 1
                        y_label = gene
                    else
                        y_label = ""
                    end
                    if j == 1
                        title_name = group
                    else
                        title_name = ""
                    end
                    ax1 = MK.Axis(fig[j,i]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
                            xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
                            xgridvisible = false, ygridvisible = false,yreversed=false, title = title_name, 
                            titlesize = 26, xlabel = "", ylabel = y_label, ylabelsize = titlesize)
                    MK.scatter!(ax1, df_plt[!, x_col], df_plt[!, y_col]; 
                            color = df_plt.plt_color, strokewidth = 0, markersize = marker_size)
                    if i == length(group_names)
                        MK.Colorbar(fig[j,length(group_names)+1], label = "", colormap = c_map, width=10, limits = (0, maximum(gene_expr)))
                    end
                end
            end
        return fig
end

function sp_dim_plot(sp::IntegratedObject, anno::Union{Symbol, String}; 
    anno_color::Union{Nothing, Dict} = nothing, x_col::String = "x", y_col::String = "y", 
    cell_order::Union{Vector{String}, Nothing}=nothing, cell_highlight::Union{String, Nothing}=nothing,
    width=600, height=500, stroke_color=:transparent, bg_color=:white, 
    marker_size=2, do_legend=true, alpha::Real = 1, nrow = 1,
    legend_size = 10, legend_fontsize = 16, legend_ncol = 1
    )
    anno_df=deepcopy(sp.spmetaData)
    anno_df[!, anno] = string.(anno_df[!, anno])
    if isa(anno, String)
        anno=Symbol(anno)
    end
    if isa(anno_color, Nothing)
        cell_anno=unique(anno_df[!,anno])
        c_map=Colors.distinguishable_colors(length(cell_anno), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        c_map = "#" .* hex.(c_map)
        anno_color=Dict(cell_anno .=> c_map)
    end
    anno_df=DataFrames.transform(anno_df, anno => ByRow(x -> anno_color[x]) => :new_color)
    anno_df.new_color = [(i, alpha) for i in anno_df.new_color]
    group_arr = string.(anno_df[!, "dataset"])
    group_names = unique(group_arr)
    n = length(group_names)
    ncol = ceil(Int, n / nrow)
    nlast = n - (nrow - 1) * ncol 
    lgd_col = (nlast == 0) ? ncol + 1 : nlast + 1
    fig = MK.Figure(size = (width * ncol, height * nrow))
    for idx in 1:n
        group = group_names[idx]
        j = ceil(Int, idx / ncol)
        i = idx % ncol == 0 ? ncol : idx % ncol
            anno_df1 = filter(:dataset => ==(group), anno_df)
            ax1 = MK.Axis(fig[j,i]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false)
            ax2 = MK.Axis(fig[j,i]; backgroundcolor = bg_color, xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false)
            if !isa(cell_highlight, Nothing)
                anno_df2 = deepcopy(anno_df1)
                anno_df3=filter(anno => ==(cell_highlight), anno_df2)
                anno_df4=filter(anno => !=(cell_highlight), anno_df2)
                colors = unique(anno_df3.new_color)
                if idx != n
                    MK.scatter!(ax1, anno_df4[!, x_col] , anno_df4[!, y_col]; strokecolor=stroke_color, 
                        color="gray90", strokewidth=0, markersize=marker_size)
                    MK.scatter!(ax1, anno_df3[!, x_col] , anno_df3[!, y_col]; strokecolor=stroke_color, 
                        color=colors[1], strokewidth=0, markersize=marker_size)
                else
                    if do_legend
                        MK.scatter!(ax2, anno_df3[!, x_col] , anno_df3[!, y_col]; strokecolor=stroke_color, visible=false,
                            color=colors[1], strokewidth=0, markersize=2*legend_size, label=cell_highlight)
                        MK.scatter!(ax1, anno_df4[!, x_col] , anno_df4[!, y_col]; strokecolor=stroke_color, 
                            color="gray90", strokewidth=0, markersize=marker_size, label=cell_highlight)
                        MK.scatter!(ax1, anno_df3[!, x_col] , anno_df3[!, y_col]; strokecolor=stroke_color, 
                            color=colors[1], strokewidth=0, markersize=marker_size, label=cell_highlight)
                    else
                        MK.scatter!(ax1, anno_df4[!, x_col] , anno_df4[!, y_col]; strokecolor=stroke_color, 
                            color="gray90", strokewidth=0, markersize=marker_size)
                        MK.scatter!(ax1, anno_df3[!, x_col] , anno_df3[!, y_col]; strokecolor=stroke_color, 
                            color=colors[1], strokewidth=0, markersize=marker_size)
                    end
                    if do_legend
                        MK.Legend(fig[nrow, lgd_col], ax2, framecolor=:white, labelsize=legend_fontsize, nbanks=legend_ncol)
                    end
                end
            else
                if isa(cell_order, Nothing)
                    cell_anno=unique(anno_df1[!,anno])
                else
                    cell_anno=cell_order
                end
                if idx != n
                    for i in cell_anno
                        anno_df2=filter(anno => ==(i), anno_df1)
                        x_ax = anno_df2[!, x_col]
                        y_ax = anno_df2[!, y_col]
                        colors = unique(anno_df2.new_color)
                        MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                                color=colors[1], strokewidth=0, markersize=marker_size)
                    end
                else
                    for i in cell_anno
                        anno_df2=filter(anno => ==(i), anno_df1)
                        x_ax = anno_df2[!, x_col]
                        y_ax = anno_df2[!, y_col]
                        colors = unique(anno_df2.new_color)
                        if do_legend
                            MK.scatter!(ax2, x_ax , y_ax; strokecolor=stroke_color, visible=false,
                                color=colors[1], strokewidth=0, markersize=2*legend_size, label=i)
                            MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                                color=colors[1], strokewidth=0, markersize=marker_size, label=i)
                        else
                            MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                                color=colors[1], strokewidth=0, markersize=marker_size)
                        end
                    end
                    if do_legend
                        MK.Legend(fig[nrow, lgd_col], ax2, framecolor=:white, labelsize=legend_fontsize, nbanks=legend_ncol)
                    end
                end
            end
        
    end
    return fig
end
