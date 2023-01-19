
function sp_dot_plot(sp::Union{CartanaObject, VisiumObject}, genes::Union{Vector, String},
    cluster::Union{Symbol, String};expr_cutoff::Union{Float64, Int64}=0, split_by::Union{String, Nothing}=nothing,
    x_title="Gene",y_title="Cell type", cell_order::Union{Vector, String, Nothing}=nothing,
    fontsize::Int64=12, color_scheme::String="yelloworangered",reverse_color::Bool=false,
    fig_height::Union{String, Int64}=400, fig_width::Union{String, Int64}=400)
    if isdefined(sp, :normCount)
        norm_counts=sp.normCount
    else
        error("Please normalize the data first!")
    end
    if isa(split_by, Nothing)
        all_df=DataFrame()
        for (i, gene) in enumerate(genes)
            gene_expr = subset_count(norm_counts; genes = [gene])
            gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            df = DataFrame()
            df.gene=gene_expr
            df.celltype=string.(sp.spmetaData.cell[!, cluster])
            avg_expr=DataFrames.combine(groupby(df, :celltype), :gene => mean => :avg_exp);
            perc_expr=DataFrames.combine(groupby(df, :celltype), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
            df_plt=innerjoin(avg_expr, perc_expr, on = :celltype)
            df_plt.gene.=gene
            all_df=[all_df; df_plt]
        end
        p=all_df |> @vlplot(:circle,
            x={"gene:o", title="Gene", scale={
                    domain=genes
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            y={"celltype:o", title="Cell type",
            scale={
                    domain=cell_order
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            color={"avg_exp:q",
                    scale={scheme=color_scheme,reverse=reverse_color}},
            size={"perc_exp:q", legend={symbolFillColor="transparent"}},
            height= fig_height, width=fig_width
            )
    else
        all_df=DataFrame()
        for (i, gene) in enumerate(genes)
            gene_expr = subset_count(norm_counts; genes = [gene])
            gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            df = DataFrame()
            df.gene=gene_expr
            df.celltype=string.(sp.spmetaData.cell[!, cluster])
            df.split_by = string.(sp.spmetaData.cell[!, split_by])
            avg_expr=DataFrames.combine(groupby(df, [:celltype, :split_by]), :gene => mean => :avg_exp)
            perc_expr=DataFrames.combine(groupby(df, [:celltype,:split_by]), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
            df_plt=innerjoin(avg_expr, perc_expr, on = [:celltype,:split_by])
            df_plt.gene.=gene
            all_df=[all_df; df_plt]
        end
        p=all_df |> @vlplot(:circle,
            x={"gene:o", title="Gene", scale={
                    domain=genes
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            y={"celltype:o", title="Cell type",
            scale={
                    domain=cell_order
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            color={"avg_exp:q",
                    scale={scheme=color_scheme,reverse=reverse_color}},
            size={"perc_exp:q", legend={symbolFillColor="transparent"}},
            column={:split_by, header={labelFontSize=16, title=nothing}},
            height= fig_height, width=fig_width
            )        
    end
    return p
end

function sp_feature_plot(sp::Union{CartanaObject, VisiumObject}, gene_list::Union{String, Vector{String}, Tuple{String}}; layer::String = "cells", x_col::Union{String, Symbol}="x",
    y_col::Union{String, Symbol}="y", cell_col = "cell", x_lims=nothing, y_lims=nothing, marker_size=2, order::Bool=true, scale::Bool = false,titlesize::Int64=24, 
    height::Real = 500, width::Real = 500, combine = true, img_res::String = "low",  adjust_contrast::Real = 1.0, adjust_brightness::Real = 0.3,
    color_keys=["gray94","orange","red3"], gene_colors = nothing, alpha::Real = 1, legend_fontsize = 10, do_legend=false, legend_size = 10)
    if isa(gene_list, String)
        gene_list = [gene_list]
    end
    n_rows = Int(ceil(length(gene_list) / 3))
    if length(gene_list) < 4
        n_cols = length(gene_list)
    else
        n_cols = 3
    end
    if layer === "cells"
        if isa(sp, VisiumObject)
            coord_cell = deepcopy(sp.spmetaData)
            x_col = Symbol(x_col)
            y_col = Symbol(y_col)
            rename!(coord_cell, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])
            if img_res == "high"
                scale_factor = sp.imageData.jsonParameters["tissue_hires_scalef"]
            elseif img_res == "low"
                scale_factor = sp.imageData.jsonParameters["tissue_lowres_scalef"]
            else
                error("img_res can only be \"high\" or \"low\"!")
            end
            coord_cell[!, x_col] =  coord_cell[!, x_col] .* scale_factor
            coord_cell[!, y_col] =  coord_cell[!, y_col] .* scale_factor
        else
            coord_cell=deepcopy(sp.spmetaData.cell)
        end
        if isdefined(sp, :normCount)
            norm_counts=sp.normCount
        else
            error("Please normalize the data first!")
        end
        if isa(x_lims, Nothing)
            x_lims=(minimum(coord_cell[!, x_col])-0.05*maximum(coord_cell[!, x_col]),1.05*maximum(coord_cell[!, x_col]))
        end
        if isa(y_lims, Nothing)
            y_lims=(minimum(coord_cell[!, y_col])-0.05*maximum(coord_cell[!, y_col]),1.05*maximum(coord_cell[!, y_col]))
        end
        c_map = ColorSchemes.ColorScheme([parse(Colorant, color_keys[1]),parse(Colorant, color_keys[2]),parse(Colorant, color_keys[3])])
        fig = MK.Figure(resolution = (width * n_cols, height * n_rows))
        for (i, gene) in enumerate(gene_list)
            gene_expr = subset_count(norm_counts; genes = [gene])
            gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            if scale
                gene_expr = unit_range_scale(gene_expr)
            end
            df = DataFrame()
            df.gene_expr = gene_expr
            coord_cell[!, cell_col] = string.(coord_cell[!, cell_col])
            df[!, cell_col] = string.(coord_cell[!, cell_col])
            df_plt = innerjoin(df, coord_cell, on = cell_col)
            df_plt.gene .= gene
            if sum(gene_expr) > 0.0
                colors = get(c_map, gene_expr, :extrema)
                plt_color = "#" .* hex.(colors)
                plt_color = [(i, alpha) for i in plt_color]
                df_plt.plt_color = plt_color
                if order
                    df_plt = sort(df_plt,:gene_expr)
                end
            else
                plt_color = repeat([color_keys[1]], length(gene_expr))
                df_plt.plt_color = plt_color
            end
            n_row = Int(ceil(i/3))
            if i < 4
                n_col1 = 2i-1
                n_col2 = 2i
            else
                n_col1 = 2*(i-3*(n_rows-1))-1
                n_col2 = 2*(i-3*(n_rows-1))
            end
            ax1 = MK.Axis(fig[n_row,n_col1]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
            xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
            xgridvisible = false, ygridvisible = false,yreversed=false, title = gene_list[i], 
            titlesize = titlesize, xlabel = "", ylabel = "", 
            xlabelsize = titlesize -4, ylabelsize = titlesize -4)
            if isa(sp, VisiumObject)
                if img_res == "high"
                    img = deepcopy(sp.imageData.highresImage)
                else
                    img = deepcopy(sp.imageData.lowresImage)
                end
                img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
                MK.image!(ax1, img2)
            end
            MK.xlims!(MK.current_axis(), x_lims)
            MK.ylims!(MK.current_axis(), y_lims)
            MK.scatter!(ax1, df_plt[!, x_col], df_plt[!, y_col]; color = df_plt.plt_color, strokewidth = 0, markersize = marker_size)
            MK.Colorbar(fig[n_row,n_col2], label = "", colormap = c_map, width=10, limits = (0, maximum(gene_expr)))
        end
        MK.current_figure()
    elseif layer === "transcripts"
            if isa(sp, VisiumObject)
                error("Visium object doesn't support transcript plot.")
            end
            coord_molecules=deepcopy(sp.spmetaData.molecule)
            if isa(x_lims, Nothing)
                x_lims=(minimum(coord_molecules[!, x_col])-0.05*maximum(coord_molecules[!, x_col]),1.05*maximum(coord_molecules[!, x_col]))
            end
            if isa(y_lims, Nothing)
                y_lims=(minimum(coord_molecules[!, y_col])-0.05*maximum(coord_molecules[!, y_col]),1.05*maximum(coord_molecules[!, y_col]))
            end
            if combine
                if isa(gene_colors, Nothing)
                    c_map=Colors.distinguishable_colors(length(gene_list), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
                    c_map = "#" .* hex.(c_map)
                else
                    c_map = gene_colors
                end
                gene_color=Dict(gene_list .=> c_map)
                gene_color["others"] = color_keys[1]
                from = collect(keys(gene_color))
                to = collect(values(gene_color))
                df_plt=DataFrames.transform(coord_molecules, :gene => ByRow(name -> name ∈ gene_list ? name : "others") => :new_gene)
                df_plt = map_values(df_plt, :new_gene, :forcolor, from, to)
                df_plt.new_gene = string.(df_plt.new_gene)
                df_plt.forcolor = [(i, alpha) for i in df_plt.forcolor]
                df_plt1 = filter(:gene => x -> x ∈ gene_list, df_plt)
                df_plt2 = filter(:gene => x -> x ∉ gene_list, df_plt)
                fig = MK.Figure(resolution = (width, height))
                ax1 = MK.Axis(fig[1,1]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
                    xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
                    xgridvisible = false, ygridvisible = false,yreversed=false, 
                    titlesize = titlesize, xlabel = "", ylabel = "", 
                    xlabelsize = titlesize -4, ylabelsize = titlesize -4)
                ax2 = MK.Axis(fig[1,1]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
                    xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
                    xgridvisible = false, ygridvisible = false,yreversed=false, 
                    titlesize = titlesize, xlabel = "", ylabel = "", 
                    xlabelsize = titlesize -4, ylabelsize = titlesize -4)
                ax3 = MK.Axis(fig[1,1]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
                    xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
                    xgridvisible = false, ygridvisible = false,yreversed=false,  
                    titlesize = titlesize, xlabel = "", ylabel = "", 
                    xlabelsize = titlesize -4, ylabelsize = titlesize -4)
                all_genes = ["others"; from[from .!= "others"]]
                all_colors = [color_keys[1]; to[to .!= color_keys[1]]]
                for (gene, ann_color) in zip(all_genes, all_colors)
                        x_ax = df_plt[!, x_col][df_plt.new_gene .== gene]
                        y_ax = df_plt[!, y_col][df_plt.new_gene .== gene]
                    if do_legend
                        MK.scatter!(ax1, x_ax , y_ax; color = ann_color, strokewidth = 0, markersize = legend_size, label = gene)
                        MK.scatter!(ax2, x_ax, y_ax; color = :white, strokewidth = 0, markersize = legend_size * 2, label = gene)
                        MK.scatter!(ax3, x_ax, y_ax; color = ann_color, strokewidth = 0, markersize = marker_size, label = gene)
                        MK.Legend(fig[1, 2], ax1, framecolor=:white, labelsize=legend_fontsize)
                    else
                        MK.scatter!(ax1, x_ax , y_ax; color = ann_color, strokewidth = 0, markersize = marker_size)
                    end
                end
                MK.xlims!(ax1, x_lims)
                MK.ylims!(ax1, y_lims)
                MK.xlims!(ax2, x_lims)
                MK.ylims!(ax2, y_lims)
                MK.xlims!(ax3, x_lims)
                MK.ylims!(ax3, y_lims)
                MK.current_figure()
            else
                fig = MK.Figure(resolution = (width * n_cols, height * n_rows))
                for (i, gene) in enumerate(gene_list)
                    n_row = Int(ceil(i/3))
                    if i < 4
                        n_col = i
                    else
                        n_col = i-3*(n_rows-1)
                    end
                    df_plt = DataFrames.transform(coord_molecules, :gene => ByRow(name -> name == gene ? color_keys[3] : color_keys[1]) => :forcolor)
                    df_plt1 = filter(:forcolor => x -> x == color_keys[1], df_plt)
                    df_plt2 = filter(:forcolor => x -> x == color_keys[3], df_plt)
                    df_plt.forcolor = [(i, alpha) for i in df_plt.forcolor]
                    df_plt1.forcolor = [(i, alpha) for i in df_plt1.forcolor]
                    df_plt2.forcolor = [(i, alpha) for i in df_plt2.forcolor]
                    ax1 = MK.Axis(fig[n_row,n_col]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
                    xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
                    xgridvisible = false, ygridvisible = false,yreversed=false, title = gene_list[i], 
                    titlesize = titlesize, xlabel = "", ylabel = "", 
                    xlabelsize = titlesize -4, ylabelsize = titlesize -4)
                    if order
                        MK.scatter!(ax1, df_plt1[!, x_col], df_plt1[!, y_col]; color = df_plt1.forcolor, strokewidth = 0, markersize = marker_size)
                        MK.scatter!(ax1, df_plt2[!, x_col], df_plt2[!, y_col]; color = df_plt2.forcolor, strokewidth = 0, markersize = marker_size)
                    else
                        MK.scatter!(ax1, df_plt[!, x_col], df_plt[!, y_col]; color = df_plt.forcolor, strokewidth = 0, markersize = marker_size)
                    end
                    MK.xlims!(ax1, x_lims)
                    MK.ylims!(ax1, y_lims)
                end
                MK.current_figure()
            end
    else
        error("Layer must be \"cells\" or \"transcripts\"")
    end
end

function plot_gene_polygons(sp::CartanaObject, gene_list::Union{String, Vector{String}, Tuple{String}};
    color_keys::Union{Vector{String}, Tuple{String,String,String}}=["gray94","orange","red3"],
    x_lims=nothing, y_lims=nothing,width=900,height=1000,stroke_width=0,stroke_color="black",
    titlesize::Int64=24, scale::Bool = false
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
    norm_count=sp.polynormCount
    polygons=sp.polygonData
    if isa(x_lims, Nothing)
        x_lims=(minimum(sp.spmetaData.cell.x)-0.05*maximum(sp.spmetaData.cell.x),1.05*maximum(sp.spmetaData.cell.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(sp.spmetaData.cell.y)-0.05*maximum(sp.spmetaData.cell.y),1.05*maximum(sp.spmetaData.cell.y))
    end
    select_fov = filter([:x, :y] => (x, y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], sp.spmetaData.cell)
    subset_poly = filter(:mapped_cell => x -> x ∈ select_fov.cell, sp.spmetaData.polygon)
    polygon_num = subset_poly.polygon_number
    polygons = polygons[polygon_num]
    c_map = ColorSchemes.ColorScheme([parse(Colorant, color_keys[1]),parse(Colorant, color_keys[2]),parse(Colorant, color_keys[3])])
    fig = MK.Figure(resolution = (width * n_cols, height * n_rows))
    for (i, gene) in enumerate(gene_list)
        gene_expr = subset_count(norm_count; genes = [gene])
        gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
        if scale
            gene_expr = unit_range_scale(gene_expr)
        end
        colors = get(c_map, gene_expr, :extrema)
        plt_color="#" .* hex.(colors)
        plt_color = plt_color[polygon_num]
        n_row = Int(ceil(i/3))
        if i < 4
            n_col1 = 2i-1
            n_col2 = 2i
        else
            n_col1 = 2*(i-3*(n_rows-1))-1
            n_col2 = 2*(i-3*(n_rows-1))
        end
        ax1 = MK.Axis(fig[n_row,n_col1]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
        xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false,yreversed=false, title = gene_list[i], 
        titlesize = titlesize, xlabel = "", ylabel = "", 
        xlabelsize = titlesize -4, ylabelsize = titlesize -4)
        MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in polygons]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width,label="")
        MK.xlims!(MK.current_axis(), x_lims)
        MK.ylims!(MK.current_axis(), y_lims)
        MK.Colorbar(fig[n_row,n_col2], label = "", colormap = c_map, width=10, limits = (0, maximum(gene_expr)))
    end
        MK.current_figure()
end

function plot_cell_polygons(sp::CartanaObject, anno;
    cell_types::Union{String, Int64, Vector, Tuple, Nothing}=nothing, cell_colors::Union{Nothing, String, Vector, Tuple} = nothing,
    bg_color = "gray90", x_lims=nothing, y_lims=nothing,width = 500, height = 500,stroke_width=0, stroke_color="black",
    legend_fontsize = 10, do_legend=false, legend_size = 10 
    )
    if isa(x_lims, Nothing)
        x_lims=(minimum(sp.spmetaData.cell.x)-0.05*maximum(sp.spmetaData.cell.x),1.05*maximum(sp.spmetaData.cell.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(sp.spmetaData.cell.y)-0.05*maximum(sp.spmetaData.cell.y),1.05*maximum(sp.spmetaData.cell.y))
    end
    if isa(cell_types, String)
        cell_types = [cell_types]
    end
    if isa(cell_colors, String)
        cell_colors = [cell_colors]
    end
    anno_df=deepcopy(sp.spmetaData.polygon)
    polygons=deepcopy(sp.polygonData)
    if isa(anno, String)
        anno=Symbol(anno)
    end
    all_celltypes = unique(anno_df[!,anno])
    if isa(cell_types, Nothing)
        cell_types = all_celltypes
    end
    other_cells = setdiff(all_celltypes, cell_types)
    other_color = Dict(other_cells .=> repeat([bg_color], length(other_cells)))
    if isa(cell_colors, Nothing)
        c_map= Colors.distinguishable_colors(length(cell_types), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        c_map = "#" .* hex.(c_map)
        cell_color=Dict(cell_types .=> c_map)
        anno_color = merge(cell_color, other_color)
    else
        if length(cell_types) !== length(cell_colors)
            error("The number of colors must equal to the number of cell types!")
        end
        cell_color=Dict(cell_types .=> cell_colors)
        anno_color = merge(cell_color, other_color)
    end
    anno_df = DataFrames.transform(anno_df, anno => ByRow(x -> anno_color[x]) => :new_color)
    plt_color = anno_df.new_color
    select_fov = filter([:x, :y] => (x, y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], sp.spmetaData.cell)
    subset_poly = filter(:mapped_cell => x -> x ∈ select_fov.cell, sp.spmetaData.polygon)
    polygon_num = subset_poly.polygon_number
    polygons = polygons[polygon_num]
    plt_color = plt_color[polygon_num]
    fig = MK.Figure(resolution=(width, height))
    ax1 = MK.Axis(fig[1,1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
        xgridvisible = false,ygridvisible = false)
        if do_legend
            cells = string.(collect(keys(cell_color)))
            colors = collect(values(cell_color))
            for (cell1, color1) in zip(cells, colors)
                select_fov1 = filter(anno => x -> x == cell1, select_fov)
                MK.scatter!(ax1, select_fov1.x , select_fov1.y; color = color1, strokewidth = 0.5,strokecolor=stroke_color, markersize = legend_size, label = cell1)
            end
            MK.Legend(fig[1, 2], ax1, framecolor=:white, labelsize=legend_fontsize)
        end
        MK.poly!(ax1, [MK.Point2.(eachrow(p)) for p in polygons]; strokecolor=stroke_color, color=plt_color, strokewidth=stroke_width)
        MK.xlims!(ax1, x_lims)
        MK.ylims!(ax1, y_lims)
    return MK.current_figure()
end

function plot_transcript_polygons(sp::CartanaObject; 
    gene_list::Union{Vector, Symbol, String, Nothing}=nothing, 
    gene_colors::Union{Vector, Symbol, String, Nothing}=nothing, 
    bg_color::Union{Vector, Symbol, String}="gray95",
    stroke_color::Union{Vector, Symbol, String}="black",
    marker_size = 2, stroke_width=0.1,
    canvas_size=(900,1000),x_lims=nothing, y_lims=nothing, 
    anno::Union{Symbol, String, Nothing}=nothing,
    ann_colors::Union{Nothing, Dict}=nothing
)
    if isa(gene_list, Nothing)
        error("Please enter a gene or a list of gene. \"gene_list\" can not be empty.")
    end
    if isa(gene_colors, Nothing)
        error("Please define the colors to label the genes. \"gene_colors\" can not be empty.")
    end
    df_spatial=sp.spmetaData.molecule
    polygons=sp.polygonData
    df_spatial[!, anno] = string.(df_spatial[!, anno])
    other_genes=unique(df_spatial.gene[Not(in.(df_spatial.gene, [Set(gene_list)]))])
    other_colors=repeat([(bg_color,0.1)],length(other_genes))
    all_genes=[gene_list; other_genes]
    all_colors=[gene_colors; other_colors]
    map_color=Dict(all_genes .=> all_colors)
    df_spatial=DataFrames.transform(df_spatial, :gene => ByRow(x -> map_color[x]) => :new_color)
    if isa(x_lims, Nothing)
        x_lims=(minimum(df_spatial.x)-0.05*maximum(df_spatial.x),1.05*maximum(df_spatial.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(df_spatial.y)-0.05*maximum(df_spatial.y),1.05*maximum(df_spatial.y))
    end 
    if isa(anno, Nothing)
        error("Please define the cell type annotation column. \"anno\" can not be empty.")
    end
    if isa(anno, String)
        anno = Symbol(anno)
    end
    anno_df=sp.spmetaData.polygon
    if isa(ann_colors, Nothing)
        cell_anno=unique(anno_df[!,anno])
        c_map=Colors.distinguishable_colors(length(cell_anno), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        c_map = "#" .* hex.(c_map)
        ann_colors=Dict(cell_anno .=> c_map)
    end
    anno_df[!, anno] = string.(anno_df[!, anno])
    anno_df=DataFrames.transform(anno_df, anno => ByRow(x -> ann_colors[x]) => :new_color)
    df_spatial1=filter(:gene => x -> x in other_genes, df_spatial)
    colors1 = df_spatial1[!,:new_color]
    df_spatial2=filter(:gene => x -> x in gene_list, df_spatial)
    colors2 = df_spatial2[!,:new_color]
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, 
        xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, 
        yticklabelsvisible=false, xgridvisible = false, ygridvisible = false)
    MK.poly!([MK.Point2.(eachrow(p)) for p in polygons]; 
        strokecolor=stroke_color, color=anno_df.new_color, strokewidth=stroke_width,label="")
    MK.scatter!(df_spatial1.x, df_spatial1.y; color=colors1,
            strokewidth=0, markersize=marker_size*0.25)
    MK.scatter!(df_spatial2.x, df_spatial2.y; color=colors2,
            strokewidth=0, markersize=marker_size)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    return MK.current_figure()
end

function sp_dim_plot(sp::Union{CartanaObject, VisiumObject}, anno::Union{Symbol, String}; 
    anno_color::Union{Nothing, Dict} = nothing, x_col::String = "x", y_col::String = "y", cell_order::Union{Vector{String}, Nothing}=nothing,
    x_lims=nothing, y_lims=nothing,canvas_size=(5000,6000),stroke_width=0.5,stroke_color=:transparent, 
        marker_size=1, label_size=50, label_color="black", label_offset=(0,0), do_label=true, do_legend=true, alpha::Real = 1,
        legend_size = 10, legend_fontsize = 16,img_res::String = "low",  adjust_contrast::Real = 1.0, adjust_brightness::Real = 0.3
    )
    if isa(sp, VisiumObject)
        anno_df = deepcopy(sp.spmetaData)
        anno_df[!, anno] = sp.metaData[!, anno]
        x_col = Symbol(x_col)
        y_col = Symbol(y_col)
        rename!(anno_df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])
        if img_res == "high"
            scale_factor = sp.imageData.jsonParameters["tissue_hires_scalef"]
        elseif img_res == "low"
            scale_factor = sp.imageData.jsonParameters["tissue_lowres_scalef"]
        else
            error("img_res can only be \"high\" or \"low\"!")
        end
        anno_df[!, x_col] =  anno_df[!, x_col] .* scale_factor
        anno_df[!, y_col] =  anno_df[!, y_col] .* scale_factor
    else
        anno_df=deepcopy(sp.spmetaData.cell)
        anno_df[!, anno] = string.(anno_df[!, anno])
    end
    if isa(x_lims, Nothing)
        x_lims=(minimum(anno_df[!,x_col])-0.05*maximum(anno_df[!,x_col]),1.05*maximum(anno_df[!,x_col]))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(anno_df[!,y_col])-0.05*maximum(anno_df[!,y_col]),1.05*maximum(anno_df[!,y_col]))
    end
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
    fig = MK.Figure(resolution=canvas_size)
    ax1 = MK.Axis(fig[1,1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
        xgridvisible = false,ygridvisible = false);
    ax2 = MK.Axis(fig[1,1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
        xgridvisible = false,ygridvisible = false);
    ax3 = MK.Axis(fig[1,1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
        xgridvisible = false,ygridvisible = false, framecolor=:white)
    if isa(cell_order, Nothing)
        cell_anno=unique(anno_df[!,anno])
    else
        cell_anno=cell_order
    end
    if isa(sp, VisiumObject)
        if img_res == "high"
            img = deepcopy(sp.imageData.highresImage)
        else
            img = deepcopy(sp.imageData.lowresImage)
        end
        img2 = augment(img, ColorJitter(adjust_contrast, adjust_brightness))
        MK.image!(ax1, img2)
    end
    for i in cell_anno
        anno_df2=filter(anno => x-> x == i, anno_df)
        x_ax = anno_df2[!, x_col]
        y_ax = anno_df2[!, y_col]
        colors = unique(anno_df2.new_color)
        if do_legend
            if isa(sp, VisiumObject)
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=colors[1], strokewidth=0, 
                    markersize=marker_size, label=i) 
            else
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=colors[1], strokewidth=0, markersize=legend_size, label=i)
                MK.scatter!(ax2, x_ax , y_ax; strokecolor=stroke_color, 
                    color=:white, strokewidth=0, markersize=2*legend_size, label=i)
                MK.scatter!(ax3, x_ax , y_ax; strokecolor=stroke_color, 
                    color=colors[1], strokewidth=0, markersize=marker_size, label=i)
            end
            MK.Legend(fig[1, 2], ax1, framecolor=:white, labelsize=legend_fontsize)
        else
            MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                color=colors[1], strokewidth=0, markersize=marker_size)
        end
    end
    if do_label
        for i in cell_anno
            anno_df = filter([x_col, y_col] => (x,y) -> x_lims[1] < x < x_lims[2] && y_lims[1] < y < y_lims[2], anno_df)
            anno_df2 = filter(anno => x-> x == i, anno_df)
            x_ax = anno_df2[!, x_col]
            y_ax = anno_df2[!, y_col]
            MK.text!(i, position = (mean(x_ax) - label_offset[1], mean(y_ax) - label_offset[2]),align = (:center, :center),font = "Noto Sans Regular",textsize = label_size,color = label_color)
        end
    end
    MK.xlims!(ax1, x_lims)
    MK.ylims!(ax1, y_lims)
    MK.xlims!(ax2, x_lims)
    MK.ylims!(ax2, y_lims)
    MK.xlims!(ax3, x_lims)
    MK.ylims!(ax3, y_lims)
    MK.ylims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    return MK.current_figure()
end

function sp_highlight_cells(sp::Union{CartanaObject, VisiumObject}, cell_hightlight::String, anno::Union{String,Symbol};
    canvas_size=(900,1000),stroke_width::Float64=0.1, stroke_color="black", cell_color::String="red",
    marker_size=2,x_lims=nothing, y_lims=nothing)
    coord_cells=deepcopy(sp.spmetaData.cell)
    if isa(x_lims, Nothing)
        x_lims=(minimum(sp.spmetaData.cell.x)-0.05*maximum(sp.spmetaData.cell.x),1.05*maximum(sp.spmetaData.cell.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(sp.spmetaData.cell.y)-0.05*maximum(sp.spmetaData.cell.y),1.05*maximum(sp.spmetaData.cell.y))
    end
    coord_cells[!,anno]=string.(coord_cells[!,anno])
    coord_cells=DataFrames.transform(coord_cells, anno => ByRow(name -> name == cell_hightlight ? name : "others") => :newcell)
    coord_cells=DataFrames.transform(coord_cells, :newcell => ByRow(name -> name =="others" ? "gray90" : cell_color) => :newcolor)
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, 
                xticksvisible=false, xticklabelsvisible=false, 
                yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false)
    MK.scatter!(coord_cells.x, coord_cells.y; color=coord_cells.newcolor, strokewidth=stroke_width, markersize=marker_size,strokecolor=stroke_color)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
end

function sp_gene_rank(sp::CartanaObject, celltype::String, cluster::String; num_gene::Int64=20)
    gene_list=unique(sp.spmetaData.molecule.gene)
    all_df=DataFrame()
    for (i, gene) in enumerate(gene_list)
        norm_counts = sp.normCount
        gene_expr = subset_count(norm_counts; genes = [gene])
        gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
        df = DataFrame()
        df.gene=gene_expr
        df.celltype=string.(sp.spmetaData.cell[!, cluster])
        avg_expr=combine(groupby(df, :celltype), :gene => mean => :avg_exp);
        avg_expr.gene .= gene
        all_df=[all_df; avg_expr]
    end
    gene_mean=DataFrame(gene=sp.rawCount.gene_name,all_mean=vec(mean(sp.rawCount.count_mtx, dims=2)))
    df_plot=innerjoin(all_df, gene_mean, on = :gene)
    clustern=filter(:celltype => x-> x == celltype, df_plot)
    clustern.rank=clustern.avg_exp .^ 2 ./ clustern.all_mean
    clustern=sort(clustern, :rank, rev=true)
    clustern=clustern[1:num_gene,:]
    clustern |> @vlplot(:circle, 
        x={:gene,scale={domain=clustern.gene},axis={title="Gene", grid=false}}, 
        y={:rank,axis={title="Ranking", grid=false}})
end

function sp_impute_gene_plot(sp::CartanaObject, gene::String; data_type="predicted", imp_type::String="SpaGE", c_map=nothing, x_lims=nothing,
    y_lims=nothing, canvas_size=(1000,1200), marker_size=2, order=true)
    if data_type === "predicted"
        if imp_type === "tangram"
            gene_count = sp.imputeData.tg_data
        elseif imp_type === "SpaGE"
            gene_count = sp.imputeData.spage_data
        elseif imp_type === "gimVI"
            gene_count = sp.imputeData.gimvi_data
        else
            error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
        end
    elseif data_type === "measured"
        gene_count=sp.normCount
    else
        error("data argument can only be \"predicted\" or \"measured\"")
    end
    df3=deepcopy(sp.spmetaData.cell)
    all_genes = gene_count.gene_name
    if !(gene in all_genes)
        gene_expr=zeros(ncol(gene_count.count_mtx))
    else
        gene_expr = subset_count(gene_count; genes = [gene])
        gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
        gene_expr=scale_data(gene_expr)
    end
    df3.gene = gene_expr
    if isa(x_lims, Nothing)
        x_lims=(minimum(sp.spmetaData.cell.x)-0.05*maximum(sp.spmetaData.cell.x),1.05*maximum(sp.spmetaData.cell.x))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(sp.spmetaData.cell.y)-0.05*maximum(sp.spmetaData.cell.y),1.05*maximum(sp.spmetaData.cell.y))
    end
    if isa(c_map, Nothing)
        c_map = ColorSchemes.ColorScheme([colorant"gray94", colorant"pink",colorant"red", colorant"red3"])
    end
    if sum(gene_expr) !==0.0
        colors = get.(Ref(c_map), (gene_expr .- minimum(gene_expr)) ./ maximum(gene_expr))
        plt_color="#" .* hex.(colors)
        df3.color = plt_color
    else
        df3.color .= "gray94"
    end
    if order
        sort!(df3, :gene)
    end
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
        xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
        xgridvisible = false,ygridvisible = false);
    MK.scatter!(df3.x, df3.y; strokecolor=:transparent, color=df3.color, strokewidth=0,label="", markersize=marker_size)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
end

function sp_impute_gene_plot(impute_list::Vector{CartanaObject}, gene_list::Vector{String}; data_type="predicted", imp_type::String="SpaGE",
    c_map=nothing, marker_size = 2, order=false, canvas_size=(500, 550))
    fig = MK.Figure(resolution=(canvas_size[1] * length(gene_list) ,canvas_size[2] * length(impute_list)))
    for j in 1:length(gene_list)
        gene=gene_list[j]
        all_expr=[]
        for i in 1:length(impute_list)
            if data_type === "predicted"
                if imp_type === "tangram"
                    impute_data = impute_list[i].imputeData.tg_data
                elseif imp_type === "SpaGE"
                    impute_data = impute_list[i].imputeData.spage_data
                elseif imp_type === "gimVI"
                    impute_data = impute_list[i].imputeData.gimvi_data
                else
                    error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
                end
                gene_expr= subset_count(impute_data; genes = [gene])
                gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            elseif data_type === "measured"
                gene_expr=subset_count(impute_list[i].normCount; genes = [gene])
                gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            else
                error("data_type can only be \"predicted\" or \"measured\"")
            end
            all_expr=[all_expr;gene_expr]
        end
        all_expr=Float64.(all_expr)
        all_expr=scale_data(all_expr;  perc=0.0)
        if c_map===nothing
            c_map = ColorSchemes.ColorScheme([colorant"gray94",colorant"pink", colorant"red", colorant"red3"])
        end
        colors = get.(Ref(c_map), (all_expr .- minimum(all_expr)) ./ maximum(all_expr))
        plt_color="#" .* hex.(colors)
        segments=[ncol(impute_list[i].normCount.count_mtx)-1 for i in 1:length(impute_list)]
        segments=cumsum(segments)
        seg_all=[]
        for i in 1:length(segments)
            if i == 1
            seg_all=[seg_all; [1:segments[i]]]
            else
            seg_all=[seg_all; [(segments[i-1]+1):(segments[i])]]
            end
        end
        for i in 1:length(impute_list)
            data_plt = deepcopy(impute_list[i].spmetaData.cell)
            data_plt.color = plt_color[seg_all[i]]
            data_plt.gene_expr = all_expr[seg_all[i]]
            if order
                sort!(data_plt, :gene_expr)
            end
            ax = MK.Axis(fig[i, j]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false,
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false)
            MK.scatter!(ax, data_plt.x, data_plt.y; 
                color=data_plt.color, strokewidth=0,markersize=marker_size)
        end
    end
    MK.current_figure()
end

function plot_fov(sp::CartanaObject, n_fields_x::Int64, n_fields_y::Int64; 
    x_col::Union{String, Symbol}="x", y_col::Union{String, Symbol}="y", group_label::Union{Nothing, String}=nothing, 
    canvas_size=(4000,4000), cell_highlight::Union{Nothing, String, Number}=nothing, shield::Bool= false, marker_size::Union{Int64, Float64}=2)
    df = sp.spmetaData.cell
    pts, centroids=split_field(df, n_fields_x, n_fields_y)
    centroids=convert.(Tuple{Float64, Float64},centroids)
    x_lims=(minimum(df[!, x_col])-0.05*maximum(df[!, x_col]),1.05*maximum(df[!, x_col]))
    y_lims=(minimum(df[!, y_col])-0.05*maximum(df[!, y_col]),1.05*maximum(df[!, y_col]))
    label_size= 800 * 16/length(centroids)
    if label_size < 35
       label_size = 35
    elseif label_size > 1600
        label_size = 1600
    end
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false);
    if isa(group_label, Nothing) && isa(cell_highlight, Nothing)
        MK.scatter!(df[!,x_col],df[!, y_col]; strokecolor="black", color=:gray98, strokewidth=0.5,label="", markersize=marker_size)
    elseif isa(group_label, Nothing) && isa(cell_highlight, String)
        error("Please indicate the group name that contains cell type info!")
    elseif isa(group_label, String) && isa(cell_highlight, Nothing)
        error("Please indicate the cell type name to be highlighted!")
    else
        df[!,group_label]=string.(df[!,group_label])
        df=DataFrames.transform(df, group_label => ByRow(name -> name == cell_highlight ? name : "others") => :newcell)
        df=DataFrames.transform(df, :newcell => ByRow(name -> name =="others" ? "gray98" : "black") => :newcolor)
        MK.scatter!(df[!, x_col],df[!, y_col]; strokecolor="black", color=df.newcolor, strokewidth=0.5,label="", markersize=marker_size)
    end
    if shield
        label_color=:yellow1
        bg_color=(:blue, 0.3)
        font_style="Noto Sans Bold"
    else
        label_color= cgrad(:darkrainbow)[LinRange(0, 1, length(centroids))]
        bg_color=:transparent
        font_style="Noto Sans Regular"
    end
    MK.poly!([p for p in pts]; color = bg_color, strokecolor = :black, strokewidth = 3)
    MK.text!(string.(1:length(centroids)),position = centroids,align = (:center, :center),font = font_style,textsize = label_size,color = label_color)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
end

function plot_point(sp::Union{CartanaObject, VisiumObject}, pt::Vector{Float64}; 
    canvas_size=(4000,4000),marker_size=60, text_size=100, 
    pt_color="red", text_color="blue", label="point")
    df = sp.spmetaData.cell
    pt2=MK.Point2f0(pt[1], pt[2])
    x_lims=(minimum(df.x)-0.05*maximum(df.x),1.05*maximum(df.x))
    y_lims=(minimum(df.y)-0.05*maximum(df.y),1.05*maximum(df.y))
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false);
    MK.scatter!(df.x,df.y; strokecolor="black", color=:gray98, strokewidth=0.5,label="")
    MK.scatter!([pt[1]],[pt[2]]; strokecolor="black", color=pt_color, strokewidth=0.5,label="",markersize=marker_size)    
    MK.text!(label,position = pt2,align = (:center, :bottom),textsize = text_size,color = text_color)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
end

function plot_depth(sp::Union{CartanaObject, VisiumObject}; celltype::Union{String, Symbol} = :celltype,
    cmap=nothing, cell_select=nothing, fontsize=16, scale=0.8, markers=nothing)
        cells=sp.spmetaData.cell
        celltypes=cell_select
        if isa(celltypes, Nothing)
            cell_order=combine(groupby(cells, celltype),:depth=>mean=>:mean)
            sort!(cell_order, :mean)
            celltypes=cell_order[!, celltype]
            celltypes=reverse(celltypes)
        end
        if isa(cmap, Nothing)
            cell_colors = [cgrad(:thermal, [0.0, 1.0])[z] for z in cells.depth]
        else
            cell_colors = [cgrad(cmap, [0.0, 1.0])[z] for z in cells.depth]
        end
        fig = MK.Figure(resolution=(1200,600))
        ax1=MK.Axis(fig[1, 1]; xticklabelsize=(fontsize-4), yticklabelsize=fontsize, xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,xgridvisible = false,ygridvisible = false,title = "Cells colored by kidney depth",titlesize = fontsize)
        ax2=MK.Axis(fig[1,2]; xticklabelsize=(fontsize-4) ,yticklabelsize=fontsize, xticksvisible=true, xticklabelsvisible=true, yticksvisible=true, yticklabelsvisible=true,xgridvisible = false,ygridvisible = false, title = "Cell distribution from cortex to papilla",titlesize = fontsize,yticks = ((1:length(celltypes)) ./ scale,  celltypes))
        MK.scatter!(ax1, cells.x, cells.y; color = cell_colors, markersize = 2)
        for i in length(celltypes):-1:1
            cell_density=filter(:celltype => x -> x == celltypes[i], cells)
            cell_density=float.(cell_density.depth)
            d = MK.density!(ax2, cell_density,npoints=200, offset = i / scale,
                color = :x, colormap = cell_colors, colorrange = (0, 1),
                strokewidth = 1, strokecolor = :black)
            MK.translate!(d, 0, 0, -0.05i)
        end
        if markers !== nothing
        molecules=sp.spmetaData.molecule
        cell2=sp.spmetaData.cell.cell
        molecules=filter(:cell=> x-> x in cell2, molecules)
        from=cell2
        to=cells.depth
        molecules2=map_values(molecules, :cell, :depth, from, to)
        markers=reverse(markers)
        ax3=MK.Axis(fig[1, 3]; xticklabelsize=(fontsize-4) ,yticklabelsize=fontsize, 
            xticksvisible=true, xticklabelsvisible=true, yticksvisible=true, 
            yticklabelsvisible=true,xgridvisible = false,ygridvisible = false, 
            title = "Transcript distribution from cortex to papilla",
            titlesize = fontsize,yticks = ((1:length(markers)) ./ scale,  markers))
            for j in length(markers):-1:1
                cell_density2=filter(:gene => x -> x == markers[j], molecules2)
                cell_density2=float.(cell_density2.depth)
                f = MK.density!(ax3, cell_density2,npoints=200, offset = j / scale,
                        color = :x, colormap = cell_colors, colorrange = (0, 1),
                        strokewidth = 1, strokecolor = :black)
                MK.translate!(f, 0, 0, -0.05j)
            end
        end
        MK.current_figure()
end

function plot_depth_animation(sp::Union{CartanaObject, VisiumObject}, celltypes::Vector{String}, markers::Vector{String}; 
    group_label="celltype",gene_label="gene", cmap=nothing, bg_color="gray94",fontsize=16, scale=0.8, canvas_size=(1800,600), file_name="animation.gif", framerate=30,
    titles=["Cells colored by kidney depth","Cell distribution from cortex to papilla","Transcript distribution from cortex to papilla"])
    cells=sp.spmetaData.cell
    cells=filter(group_label => x-> x in celltypes, cells)
    molecules=sp.spmetaData.molecule
    molecules=filter(gene_label => x-> x in markers, cells)
    fig = MK.Figure(resolution=canvas_size)
    ax1=MK.Axis(fig[1, 1]; xticklabelsize=(fontsize-4), yticklabelsize=fontsize, xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,xgridvisible = false,ygridvisible = false,title = titles[1],titlesize = fontsize)
    ax2=MK.Axis(fig[1, 3]; xticklabelsize=(fontsize-4) ,yticklabelsize=fontsize, xticksvisible=true, xticklabelsvisible=true, yticksvisible=true, yticklabelsvisible=true,xgridvisible = false,ygridvisible = false, title = titles[2],titlesize = fontsize,yticks = ((1:length(celltypes)) ./ scale,  celltypes))
    ax3=MK.Axis(fig[1, 2]; xticklabelsize=(fontsize-4) ,yticklabelsize=fontsize, xticksvisible=true, xticklabelsvisible=true, yticksvisible=true, yticklabelsvisible=true,xgridvisible = false,ygridvisible = false, title = titles[3],titlesize = fontsize,yticks = ((1:length(markers)) ./ scale,  markers))
    MK.scatter!(ax1, cells.x, cells.y; color = bg_color, markersize = 2)
    celltypes=reverse(celltypes)
    for i in length(celltypes):-1:1
        celltype=celltypes[i]
        cell_density0=filter(:celltype => x -> x == celltype, cells)
        cell_density0=float.(cell_density0.depth)
        d0 = MK.density!(ax2, cell_density0,npoints=200, offset = i / scale, 
            color = bg_color, strokewidth = 1, strokecolor = :black)
        MK.translate!(d0, 0, 0, -0.05i)
    end
    markers = reverse(markers)
    for j in length(markers):-1:1
        mol_density0=filter(gene_label => x -> x == markers[j], molecules)
        mol_density0=float.(mol_density0.depth)
        f0 = MK.density!(ax3, mol_density0,npoints=200, offset = j / scale,
                color = bg_color, strokewidth = 1, strokecolor = :black)
        MK.translate!(f0, 0, 0, -0.05j)
    end
    timestamps = range(0.01, 1, length=framerate)
    MK.record(fig, file_name, timestamps; framerate = framerate) do t
            cells_t=filter(:depth => x -> x < t, cells)
            cell_stats_t=countmap(cells_t.celltype)
            depth_val = cells_t.depth
            depth_val=float.(depth_val)
            cs = ColorSchemes.ColorScheme([colorant"#0d2232",colorant"#1b326e",colorant"#473899",colorant"#6e4b8e",colorant"#935c87",colorant"#a6606c",colorant"#dc7f60",colorant"#eea152",colorant"#f2cd5d", colorant"#ecf975"])
            colors = get.(Ref(cs), (depth_val .- minimum(depth_val)) ./ maximum(depth_val))
            cell_colors="#" .* hex.(colors)
            MK.scatter!(ax1, cells_t.x, cells_t.y; color = cell_colors, markersize = 2)
            for i in 13:-1:1
                cell_density=filter(:celltype => x -> x == celltypes[i], cells)
                cell_density=float.(cell_density.depth)
                cs = ColorSchemes.ColorScheme([colorant"#0d2232",colorant"#1b326e",colorant"#473899",colorant"#6e4b8e",colorant"#935c87",colorant"#a6606c",colorant"#dc7f60",colorant"#eea152",colorant"#f2cd5d", colorant"#ecf975", "colorant" * "\"" * bg_color * "\""])
                sort!(depth_val)
                colors = get.(Ref(cs), (depth_val .- minimum(depth_val)) ./ maximum(depth_val))
                cell_colors="#" .* hex.(colors)
                d = MK.density!(ax2, cell_density,npoints=200, offset = i / 0.5, 
                    color = :x, colormap = cell_colors, colorrange = (0, t),
                    strokewidth = 1, strokecolor = :black)
                MK.translate!(d, 0, 0, -0.05i)
            end
            molecule_t=filter(:depth=> x-> x < t, molecules2)
            for j in 13:-1:1
                mol_density=filter(gene_label => x -> x == markers[j], molecules)
                mol_density=float.(mol_density.depth)
                cs = ColorSchemes.ColorScheme([colorant"#0d2232",colorant"#1b326e",colorant"#473899",colorant"#6e4b8e",colorant"#935c87",colorant"#a6606c",colorant"#dc7f60",colorant"#eea152",colorant"#f2cd5d", colorant"#ecf975", "colorant" * "\"" * bg_color * "\""])
                depth_val=molecule_t.depth
                sort!(depth_val)
                colors = get.(Ref(cs), (depth_val .- minimum(depth_val)) ./ maximum(depth_val))
                mol_colors="#" .* hex.(colors)
                f = MK.density!(ax3, mol_density,npoints=200, offset = j / 0.5,
                        color = :x, colormap = mol_colors, colorrange = (0, t),
                        strokewidth = 1, strokecolor = :black)
                MK.translate!(f, 0, 0, -0.05j)
            end
        end
end

function plot_gene_depth(sp::Union{CartanaObject, VisiumObject}, gene::String;
    c_map::Union{String, Symbol, Nothing}=nothing, cell_col="cell",
    canvas_size =(1200,300),marker_size=4,
    stroke_width=0.5,stroke_color="gray94",
    expr_cutoff=0.25,n_bins=50
)
    coord_cell=sp.spmetaData.cell
    norm_counts=sp.normCount
    gene_expr = subset_count(norm_counts; genes = [gene])
    gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
    df = DataFrame()
    df.gene_expr=gene_expr
    coord_cell[!, cell_col]=string.(coord_cell[!, cell_col])
    df.cell=string.(coord_cell[!, cell_col])
    df_plt=innerjoin(df, coord_cell, on = cell_col)
    df_plt.gene.=gene
    df_plt.depth=float.(df_plt.depth)
    df_plt.gene_expr=float.(df_plt.gene_expr)
    if isa(c_map, Nothing)
        c_map=:gist_heat
    end
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=16, yticklabelsize=16, 
                xticksvisible=false, xticklabelsvisible=false, 
                yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false, 
                ylabel = "Expression", xlabel = "Kidney depth")
    df_plt=filter(:gene_expr=>x-> x > expr_cutoff, df_plt)
    x_hist = fit(Histogram, tuple(df_plt.depth, df_plt.gene_expr), nbins=n_bins) 
    MK.contour!(x_hist,levels = 3,fillrange = true, colormap=c_map)
    MK.ylims!(MK.current_axis(), (0.25, 0.75))
    MK.current_figure()
end

function PlotInteractive(sp::Union{CartanaObject, VisiumObject}; layer::String = "cells", marker_color::Union{Symbol, String}="black", marker_size=3, plot_mode="markers")
    if layer === "cells"
        cells=sp.spmetaData.cell
        plyjs.plot(plyjs.scatter(x=cells.x, y=cells.y, mode=plot_mode, marker=plyjs.attr(size=marker_size, color=marker_color)))
    elseif layer === "molecules"
        molecules = sp.spmetaData.molecule
        plyjs.plot(plyjs.scatter(x=molecules.x, y=molecules.y, mode=plot_mode, marker=plyjs.attr(size=marker_size, color=marker_color)))
    else
        println("Layer must be \"cells\" or \"molecules\"")
    end
end

function plot_transcript_dapi(sp::CartanaObject, fov::Int64, n_fields_x::Int64, 
    n_fields_y::Int64; noise_ann=nothing,annotation=:cell,
    is_noise=nothing, draw_poly=false, marker_size=3)
    selected_view = subset_fov(sp, fov, n_fields_x, n_fields_y)
    xmin = trunc(Int64,minimum(selected_view.x))
    xmax = trunc(Int64,maximum(selected_view.x))
    ymin = trunc(Int64,minimum(selected_view.y))
    ymax = trunc(Int64, maximum(selected_view.y))
    polygons = sp.polygonData
    df_spatial = sp.spmetaData.molecule
    img = sp.imageData
    df_spatial = filter([:x, :y]=> (x,y) -> xmin < x < xmax && ymin < y < ymax, df_spatial);
    df_spatial = filter(:is_noise=> x -> x ==0, df_spatial)
    poly = sp.imageData.polygon
    cells = df_spatial.cell
    poly.mapped_cell = Int.(poly.mapped_cell)
    poly = filter(:mapped_cell=> x-> x in cells, poly)
    polygons = polygons[poly.polygon_number]
    img2 = img[ymin:ymax,xmin:xmax]'
    plt_x = df_spatial.x .- xmin
    plt_y = df_spatial.y .- ymin
    fig = MK.Figure(resolution=(500,500))
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, 
        xticksvisible = false, xticklabelsvisible=false, backgroundcolor = :black,
        yticksvisible = false, yticklabelsvisible=false, xgridvisible = false,ygridvisible = false );
    annotation = df_spatial[!,annotation]
    ann_vals = annotation[annotation .!= noise_ann] |> unique |> sort
    c_map = Colors.distinguishable_colors(length(ann_vals), 
        Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
    MK.image!(img2)
    for (color, ann) in zip(c_map, ann_vals)
            MK.scatter!(df_spatial.x[annotation .== ann] .- xmin, df_spatial.y[annotation .== ann] .- ymin;
            strokewidth=0, markersize=marker_size, label=ann, color=(color,1))
    end
    if draw_poly
        MK.poly!([MK.Point2.(eachrow(p .- [xmin ymin])) for p in polygons]; strokecolor=("cyan",0.7), color="transparent", strokewidth=2,label="")
    end
    MK.xlims!(MK.current_axis(), (minimum(plt_x), maximum(plt_x)))
    MK.ylims!(MK.current_axis(), (minimum(plt_y), maximum(plt_y)))
    MK.current_figure()
end

function compare_gene_imputation(sp1::CartanaObject,sp2::CartanaObject, gene_list::Union{Vector, String},
    cluster::Union{Symbol, String}; sp1_name::String ="sp1", sp2_name::String="sp2",
    assay_use::String="measured",expr_cutoff::Union{Float64, Int64}=0, legend_min::Union{Float64, Int64}=0, legend_max::Union{Float64, Int64}=1, 
    x_title="Gene",y_title="Cell type", cell_order::Union{Vector, String, Nothing}=nothing,
    fontsize::Int64=12, color_range::Vector=["white", "ivory","gold","orange","tomato","red"],
    fig_height::Union{String, Int64}=400, fig_width::Union{String, Int64}=400)
    all_df=DataFrame()
    if assay_use === "measured"
        ct_mtx1 = deepcopy(sp1.normCount)
        ct_mtx2 = deepcopy(sp2.normCount)
    elseif assay_use === "predicted"
            if imp_type === "tangram"
                ct_mtx1 = sp1.imputeData.tg_data
                ct_mtx2 = sp2.imputeData.tg_data
            elseif imp_type === "SpaGE"
                ct_mtx1 = sp1.imputeData.spage_data
                ct_mtx2 = sp2.imputeData.spage_data
            elseif imp_type === "gimVI"
                ct_mtx1 = sp1.imputeData.gimvi_data
                ct_mtx2 = sp2.imputeData.gimvi_data
            else
                error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
            end
    else
        error("assay_use can only be \"measured\" or \"predicted\"")
    end
        all_df=DataFrame()
        for (i, gene) in enumerate(gene_list)
            gene_expr= subset_count(ct_mtx1; genes = [gene])
            gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            df = DataFrame()
            df.gene=gene_expr
            df.celltype=string.(sp1.spmetaData.cell[!, cluster])
            avg_expr1=combine(groupby(df, :celltype), :gene => mean => :avg_exp)
            avg_expr1.group .= sp1_name
            gene_expr= subset_count(ct_mtx2; genes = [gene])
            gene_expr = (vec ∘ collect)(Float64.(gene_expr.count_mtx))
            df = DataFrame()
            df.gene=gene_expr
            df.celltype=string.(sp2.spmetaData.cell[!, cluster])
            avg_expr2=combine(groupby(df, :celltype), :gene => mean => :avg_exp)
            avg_expr2.group .= sp2_name
            avg_expr = [avg_expr1; avg_expr2]
            if scale
                avg_expr.avg_exp= unit_range_scale(avg_expr.avg_exp)
            end
            perc_expr=combine(groupby(df, :celltype), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
            df_plt=innerjoin(avg_expr, perc_expr, on = :celltype)
            df_plt.gene.=gene
            all_df=[all_df; df_plt]
        end
        p=all_df |> @vlplot(:rect,
            y={"gene:o", title="Gene", scale={
                    domain=gene_list
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            x={"celltype:o", title="Cell type",
               scale={
                    domain=cell_order
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            color={"avg_exp:q",
                    scale={domainMin=legend_min, domainMax=legend_max, range=color_range}},
            column={:group, header={labelFontSize=16, title=nothing}},
            height= fig_height, width=fig_width
            )
        return p
end

function plot_heatmap(sp::CartanaObject, gene_list::Union{Vector, String},
    cluster::Union{Symbol, String};assay_use::String="measured",expr_cutoff::Union{Float64, Int64}=0, split_by::Union{String, Nothing}=nothing,
    x_title="Gene",y_title="Cell type", cell_order::Union{Vector, String, Nothing}=nothing,
    fontsize::Int64=12, color_scheme::String="yelloworangered",reverse_color::Bool=false,scale::Bool=false,
    fig_height::Union{String, Int64}=400, fig_width::Union{String, Int64}=400)
    all_df=DataFrame()
    if assay_use === "measured"
        ct_mtx = deepcopy(sp.normCount)
    elseif assay_use === "predicted"
        if imp_type === "tangram"
            ct_mtx = sp.imputeData.tg_data
        elseif imp_type === "SpaGE"
            ct_mtx = sp.imputeData.spage_data
        elseif imp_type === "gimVI"
            ct_mtx = sp.imputeData.gimvi_data
        else
            error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
        end
    else
        error("assay_use can only be \"measured\" or \"predicted\"")
    end
    if isa(split_by, Nothing)
        all_df=DataFrame()
        for (i, gene) in enumerate(gene_list)
            gene_expr= subset_count(ct_mtx; genes = [gene])
            gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            df = DataFrame()
            df.gene=gene_expr
            df.celltype=string.(sp.spmetaData.cell[!, cluster])
            avg_expr=combine(groupby(df, :celltype), :gene => mean => :avg_exp)
            if scale
                avg_expr.avg_exp= unit_range_scale(avg_expr.avg_exp)
            end
            perc_expr=combine(groupby(df, :celltype), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
            df_plt=innerjoin(avg_expr, perc_expr, on = :celltype)
            df_plt.gene.=gene
            all_df=[all_df; df_plt]
        end
        p=all_df |> @vlplot(:rect,
            y={"gene:o", title="Gene", scale={
                    domain=gene_list
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            x={"celltype:o", title="Cell type",
               scale={
                    domain=cell_order
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            color={"avg_exp:q",
                    scale={scheme=color_scheme,reverse=reverse_color}},
            height= fig_height, width=fig_width
            )
    else
        all_df=DataFrame()
        for (i, gene) in enumerate(gene_list)
            gene_expr= subset_count(ct_mtx; genes = [gene])
            gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
            df = DataFrame()
            if scale
                gene_expr =unit_range_scale(gene_expr)
            end
            df.gene = gene_expr
            df.celltype=string.(sp.spmetaData.cell[!, cluster])
            df.split_by = string.(sp.spmetaData.cell[!, split_by])
            avg_expr=combine(groupby(df, [:celltype, :split_by]), :gene => mean => :avg_exp)
            perc_expr=combine(groupby(df, [:celltype,:split_by]), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
            df_plt=innerjoin(avg_expr, perc_expr, on = [:celltype,:split_by])
            df_plt.gene.=gene
            all_df=[all_df; df_plt]
        end
        p=all_df |> @vlplot(:rect,
            y={"gene:o", title="Gene", scale={
                    domain=gene_list
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            x={"celltype:o", title="Cell type",
               scale={
                    domain=cell_order
                }, axis={labelFontSize=fontsize,titleFontSize=fontsize}},
            color={"avg_exp:q",
                    scale={scheme=color_scheme,reverse=reverse_color}},
            column={:split_by, header={labelFontSize=16, title=nothing}},
            height= fig_height, width=fig_width
            )        
    end
        return p
end

function overlay_visium_cartana(vs::VisiumObject, sp::CartanaObject; vs_x = "new_x", vs_y = "new_y", 
    sp_x = "new_x", sp_y = "new_y", vs_color=:red, sp_color=:blue, vs_markersize=7, 
    sp_markersize=2, vs_title="Visium", sp_title="Cartana")
    cartana_df = deepcopy(sp.spmetaData.cell)
    visium_df = deepcopy(vs.spmetaData.cell)
    fig = MK.Figure(resolution=(1800,500))
    ax1 = MK.Axis(fig[1, 1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                            xgridvisible = false, ygridvisible = false,yreversed=false, title = sp_title, titlesize=26)
    ax2 = MK.Axis(fig[1, 2]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                            xgridvisible = false, ygridvisible = false,yreversed=false, title = vs_title, titlesize=26)
    ax3 = MK.Axis(fig[1, 3]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                            xgridvisible = false, ygridvisible = false,yreversed=false, title = sp_title * "+" * vs_title, titlesize=26)
    MK.scatter!(ax1, unit_range_scale(cartana_df[!, sp_x]), unit_range_scale(cartana_df[!, sp_y]); color=(sp_color,0.5),strokewidth=0, markersize=sp_markersize)
    MK.scatter!(ax2, unit_range_scale(visium_df[!, vs_x]), unit_range_scale(visium_df[!, vs_y]); color=(vs_color,0.8),strokewidth=0, markersize=vs_markersize)
    MK.scatter!(ax3, unit_range_scale(cartana_df[!, sp_x]), unit_range_scale(cartana_df[!, sp_y]); color=(sp_color,0.5),strokewidth=0, markersize=sp_markersize)
    MK.scatter!(ax3, unit_range_scale(visium_df[!, vs_x]), unit_range_scale(visium_df[!, vs_y]); color=(vs_color,0.8),strokewidth=0, markersize=vs_markersize)
    MK.current_figure()
end

function overlay_visium_cartana_gene(vs::VisiumObject, sp::CartanaObject, gene; vs_x="new_x", vs_y="new_y", sp_x="new_x", sp_y="new_y",
    vs_color=:red, sp_color=:blue, vs_markersize=7, canvas_size=(1800,500),x_lims=nothing, y_lims=nothing,
    sp_markersize=2, vs_title="Visium", sp_title="Cartana", order=true, scale = true)
    vs_count=deepcopy(vs.normCount)
    sp_count=deepcopy(sp.normCount)
    gene_expr= subset_count(vs_count; genes = [gene])
    gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
    if scale
        gene_expr=unit_range_scale(gene_expr)
    end
    df1=DataFrame()
    df1.gene_expr = gene_expr
    c_map = ColorSchemes.ColorScheme([colorant"gray96",colorant"red",colorant"red3"])
    colors = get.(Ref(c_map), (gene_expr .- minimum(gene_expr)) ./ maximum(gene_expr))
    plt_color1="#" .* hex.(colors)
    df1.color1=plt_color1
    plt_color2=[(x, 0.5) for x in plt_color1]
    df1.color2=plt_color2
    df1.x = vs.spmetaData.cell[!, vs_x]
    df1.y = vs.spmetaData.cells[!, vs_y]
    gene_expr= subset_count(sp_count; genes = [gene])
    gene_expr = (vec ∘ collect)(gene_expr.count_mtx)
    if scale
        gene_expr=unit_range_scale(gene_expr)
    end
    df2=DataFrame()
    df2.gene_expr = gene_expr
    c_map = ColorSchemes.ColorScheme([colorant"gray96",colorant"blue",colorant"blue3"])
    colors = get.(Ref(c_map), (gene_expr .- minimum(gene_expr)) ./ maximum(gene_expr))
    plt_color3="#" .* hex.(colors)
    df2.color3=plt_color3
    plt_color4=[(x, 0.5) for x in plt_color3]
    df2.color4=plt_color4
    df2.x = sp.spmetaData.cell[!, sp_x]
    df2.y = sp.spmetaData.cell[!, sp_y]
    if order
        df2 = sort(df2,:gene_expr)        
    end
    if isa(x_lims, Nothing)
        x_lims=(minimum(sp.spmetaData.cell[!, sp_x])-0.1,maximum(sp.spmetaData.cell[!, sp_x])+0.1)
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(sp.spmetaData.cell[!, sp_y])-0.1,maximum(sp.spmetaData.cell[!, sp_y])+0.1)
    end
    fig = MK.Figure(resolution=canvas_size)
    ax3 = MK.Axis(fig[1,3]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                            xgridvisible = false, ygridvisible = false, title = vs_title * "+" * sp_title, titlesize=26)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    ax1 = MK.Axis(fig[1,1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false,
                            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,ylabel=gene,ylabelsize=26,
                            xgridvisible = false, ygridvisible = false, title = sp_title, titlesize=26)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    ax2 = MK.Axis(fig[1,2]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
                            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
                            xgridvisible = false, ygridvisible = false, title = vs_title, titlesize=26)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.scatter!(ax1, df2.x, df2.y; color=df2.color3, strokewidth=0, markersize=sp_markersize)
    MK.scatter!(ax2, df1.x, df1.y; color = df1.color1, strokewidth=0, markersize=vs_markersize)
    MK.scatter!(ax3, df2.x, df2.y; color=df2.color4, strokewidth=0, markersize=sp_markersize)
    MK.scatter!(ax3, df1.x, df1.y; color=df1.color2, strokewidth=0, markersize=vs_markersize)
    MK.current_figure()
end

function plot_gene_group_spatial(sp_list::Vector{CartanaObject}, n_bin, gene_list; group_names::Union{Vector, String, Nothing}=nothing,
    color_range::Vector=["white", "ivory","gold","orange","tomato","red"],legend_min::Union{Float64, Int64}=0, legend_max::Union{Float64, Int64}=1)
    n_obj = length(sp_list)
    all_genes = DataFrame()
    for i in 1:n_obj
        bin_data = bin_gene_spatial(sp_list[i], n_bin)
        plt_df = filter(:gene => x -> x in gene_list, bin_data)
        if isa(group_names, Nothing)
            plt_df.group .= "group" * string(i)
        else
            plt_df.group .= group_names[i]
        end
        all_genes = [all_genes; plt_df]
    end
    plt_df = DataFrame()
    for i in 1:length(gene_list)
        gene1 = filter(:gene => x -> x == gene_list[i], all_genes)
        gene1.avg_exp = unit_range_scale(gene1.avg_exp)
        plt_df = [plt_df; gene1]
    end
    plt_df.bin = round.(plt_df.bin; digits= 2)
    p = plt_df |> @vlplot(:rect,
        y={"group:o", title="",scale={domain=group_names}, axis={labelFontSize=16,titleFontSize=16}},
        x={"bin:o", title="C --------------------- > P", axis={labelFontSize=0,titleFontSize=16}},
        color={"avg_exp:q",scale={domainMin=legend_min, domainMax=legend_max, range=color_range}},
        column={:gene,sort=gene_list, header={labelFontSize=16, title=nothing}},
        height= 200, width=180
        )
    return p
end