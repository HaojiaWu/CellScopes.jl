function dim_plot(sc_obj::Union{scRNAObject, VisiumObject,ImagingSpatialObject, CartanaObject, XeniumObject,scATACObject, MerfishObject, SlideseqObject, seqFishObject, STARmapObject,StereoSeqObject}; anno::Union{Symbol, String}="cluster", dim_type::String="umap",
    anno_color::Union{Nothing, Dict} = nothing, cell_order::Union{Vector{String}, Nothing}=nothing,
    x_lims=nothing, y_lims=nothing,canvas_size=(600,500),stroke_width=0.5,stroke_color=:transparent, 
        marker_size=2, label_size=20, label_color="black", label_offset=(0,0), do_label=true, do_legend=true,
        legend_size = 10, legend_fontsize = 16, legend_ncol=1)   
        dim_data,x_col, y_col = get_dim_data(sc_obj.dimReduction, dim_type)
        meta_data = sc_obj.metaData
        if isa(x_lims, Nothing)
            x_lims=(minimum(dim_data[!,x_col])-0.05*maximum(dim_data[!,x_col]),1.05*maximum(dim_data[!,x_col]))
        end
        if isa(y_lims, Nothing)
            y_lims=(minimum(dim_data[!,y_col])-0.05*maximum(dim_data[!,y_col]),1.05*maximum(dim_data[!,y_col]))
        end
        if anno in names(meta_data)
            dim_data[!, anno] = meta_data[!, anno]
            if isa(anno_color, Nothing)
                cell_anno=unique(dim_data[!,anno])
                c_map=Colors.distinguishable_colors(length(cell_anno), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
                c_map = "#" .* hex.(c_map)
                anno_color=Dict(cell_anno .=> c_map)
            end
            dim_data=DataFrames.transform(dim_data, Symbol(anno) => ByRow(x -> anno_color[x]) => :new_color)
        else
            dim_data.new_color .= "springgreen3"
        end
        fig = MK.Figure(resolution=canvas_size)
        ax1 = MK.Axis(fig[1,1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false, ygridvisible = false, xlabel = names(dim_data)[1], ylabel = names(dim_data)[2]);
        ax2 = MK.Axis(fig[1,1]; xticklabelsize=12, yticklabelsize=12, xticksvisible=false, 
            xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false,
            xgridvisible = false,ygridvisible = false);
        if isa(cell_order, Nothing)
            cell_anno=unique(dim_data[!,anno])
        else
            cell_anno=cell_order
        end
        for i in cell_anno
            MK.xlims!(MK.current_axis(), x_lims)
            MK.ylims!(MK.current_axis(), y_lims)
            anno_df2=filter(anno => ==(i), dim_data)
            x_ax = anno_df2[!, x_col]
            y_ax = anno_df2[!, y_col]
            colors = unique(anno_df2.new_color)
            if do_legend
                 MK.scatter!(ax2, x_ax , y_ax; strokecolor=stroke_color, visible=false,
                    color=string.(colors[1]), strokewidth=0, markersize=2*legend_size, label=i)
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=string.(colors[1]), strokewidth=0, markersize=marker_size, label=i)
            else
                MK.scatter!(ax1, x_ax , y_ax; strokecolor=stroke_color, 
                    color=string.(colors[1]), strokewidth=0, markersize=marker_size)
            end
        end
        if do_legend
            MK.Legend(fig[1, 2], ax2, framecolor=:white, labelsize=legend_fontsize, nbanks=legend_ncol)
        end
        if do_label
            for i in cell_anno
                anno_df2=filter(anno => ==(i), dim_data)
                x_ax = anno_df2[!, x_col]
                y_ax = anno_df2[!, y_col]
                MK.text!(i, position = (mean(x_ax) - label_offset[1], mean(y_ax) - label_offset[2]),align = (:center, :center),font = "Noto Sans Regular",fontsize = label_size,color = label_color)
            end
        end
        MK.current_figure()
end

function highlight_cells(sc_obj::Union{scRNAObject, VisiumObject, ImagingSpatialObject, CartanaObject, XeniumObject, scATACObject, MerfishObject, SlideseqObject,StereoSeqObject}, cl::String; dim_type="umap", anno::Union{String,Symbol}="cluster",
    canvas_size=(600,500),stroke_width::Float64=0.1, stroke_color="black", cell_color::String="red",
    marker_size=2,x_lims=nothing, y_lims=nothing)
    dim_data, x_col, y_col = get_dim_data(sc_obj.dimReduction, dim_type)
    meta_data = sc_obj.metaData
    dim_data.cluster = meta_data.cluster
    if isa(x_lims, Nothing)
        x_lims=(minimum(dim_data[!,x_col])-0.05*maximum(dim_data[!,x_col]),1.05*maximum(dim_data[!,x_col]))
    end
    if isa(y_lims, Nothing)
        y_lims=(minimum(dim_data[!,y_col])-0.05*maximum(dim_data[!,y_col]),1.05*maximum(dim_data[!,y_col]))
    end
    dim_data = DataFrames.transform(dim_data, anno => ByRow(name -> name == cl ? name : "others") => :newcell)
    dim_data = DataFrames.transform(dim_data, :newcell => ByRow(name -> name =="others" ? "gray90" : cell_color) => :newcolor)
    fig = MK.Figure(resolution=canvas_size)
    fig[1, 1] = MK.Axis(fig; xticklabelsize=12, yticklabelsize=12, 
                xticksvisible=false, xticklabelsvisible=false, 
                yticksvisible=false, yticklabelsvisible=false,
                xgridvisible = false,ygridvisible = false, 
                xlabel = names(dim_data)[1], 
                ylabel = names(dim_data)[2])
    MK.scatter!(dim_data[!,x_col], dim_data[!,y_col]; color=dim_data.newcolor, strokewidth=stroke_width, markersize=marker_size,strokecolor=stroke_color)
    MK.xlims!(MK.current_axis(), x_lims)
    MK.ylims!(MK.current_axis(), y_lims)
    MK.current_figure()
end

function feature_plot(sc_obj::Union{scRNAObject, VisiumObject, ImagingSpatialObject, CartanaObject, XeniumObject, scATACObject, MerfishObject, SlideseqObject, seqFishObject, STARmapObject, StereoSeqObject}, genes; dim_type::String = "umap", count_type = "norm",x_lims=nothing, y_lims=nothing, marker_size=4, order=true,
    color_keys::Union{Vector{String}, Tuple{String,String,String}}=("black","yellow","red"), do_dimname::Bool=false,
        split_by::Union{String, Symbol, Nothing}=nothing, titlesize::Int64 = 24, height::Real = 500, width::Real = 500)
        dim_data, x_col, y_col = get_dim_data(sc_obj.dimReduction, dim_type)
        if count_type === "norm"
            ct_obj = subset_count(sc_obj.normCount; genes = genes)
        elseif count_type === "raw"
            ct_obj = subset_count(sc_obj.rawCount; genes = genes)
        elseif count_type === "scale"
            ct_obj = subset_count(sc_obj.scaleCount; genes = genes)
        else
            println("count_type can only be \"raw\", \"norm\" or \"scale\"!")
        end
        count_mat = ct_obj.count_mtx
        count_mat = DataFrame(count_mat, :auto)
        count_mat.gene = ct_obj.gene_name
        count_mat = permutedims(count_mat, :gene)
        count_mat.cells = sc_obj.rawCount.cell_name
        gene_data = [dim_data count_mat]
        if isa(x_lims, Nothing)
            x_lims=(minimum(gene_data[!, x_col])-0.05*maximum(gene_data[!, x_col]),1.05*maximum(gene_data[!, x_col]))
        end
        if isa(y_lims, Nothing)
            y_lims=(minimum(gene_data[!, y_col])-0.05*maximum(gene_data[!, y_col]),1.05*maximum(gene_data[!, y_col]))
        end
        c_map = ColorSchemes.ColorScheme([parse(Colorant, color_keys[1]),parse(Colorant, color_keys[2]),parse(Colorant, color_keys[3])])
        if isa(split_by, Nothing)
            n_rows = Int(ceil(length(genes) / 3))
            if length(genes) < 4
                n_cols = length(genes)
            else
                n_cols = 3
            end
            fig = MK.Figure(resolution = (width * n_cols, height * n_rows))
            for (i, gene) in enumerate(genes)
                df_plt = gene_data[!, [x_col, y_col, gene]]
                gene_expr = float.(df_plt[!, gene])
                if sum(gene_expr) !== 0.0
                    @inbounds colors = get(c_map, gene_expr, :extrema)
                    plt_color = "#" .* hex.(colors)
                    df_plt.plt_color = plt_color
                    if order
                        df_plt = sort(df_plt, Symbol(gene))
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
                    n_col1 = 2*(i-3*(n_row-1))-1
                    n_col2 = 2*(i-3*(n_row-1))
                end
                if do_dimname
                    x_label = names(df_plt)[1]
                    y_label = names(df_plt)[2]
                else
                    x_label = ""
                    y_label = ""
                end
                ax1 = MK.Axis(fig[n_row,n_col1]; xticklabelsize = 12, yticklabelsize = 12, xticksvisible = false, 
                                    xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
                                    xgridvisible = false, ygridvisible = false,yreversed=false, title = genes[i], 
                                    titlesize = titlesize, xlabel = x_label, ylabel = y_label, 
                                    xlabelsize = titlesize -4, ylabelsize = titlesize -4)
                MK.scatter!(ax1, df_plt[!, x_col], df_plt[!, y_col]; 
                            color = df_plt.plt_color, strokewidth = 0, markersize = marker_size)
                MK.Colorbar(fig[n_row,n_col2], label = "", colormap = c_map, width=10, limits = (0, maximum(gene_expr)))
            end
        else
            group_arr = string.(sc_obj.metaData[!, split_by])
            group_names = unique(group_arr)
            gene_data[!, split_by] = group_arr
            fig = MK.Figure(resolution = (width * length(group_names), height * length(genes)))
            for (i, group) in enumerate(group_names)
                for (j, gene) in enumerate(genes)
                    df_plt = gene_data[!, [x_col, y_col, gene, split_by]]
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
                    df_plt = filter(split_by => ==(group), df_plt)
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
        end
        MK.current_figure()
end

function dot_plot(sc_obj::Union{scRNAObject, VisiumObject, ImagingSpatialObject, CartanaObject, XeniumObject, scATACObject, MerfishObject, SlideseqObject, seqFishObject, STARmapObject,StereoSeqObject}, genes::Union{Vector, String},
    cluster::Union{Symbol, String};count_type = "norm" , expr_cutoff::Union{Float64, Int64}=0, split_by::Union{String, Nothing}=nothing,
    x_title="Gene", y_title = "Cell type", cell_order::Union{Vector, String, Nothing}=nothing,
    fontsize::Int64 = 12, color_scheme::String="yelloworangered",reverse_color::Bool=false,
    height::Union{String, Int64}=400, width::Union{String, Int64}=400)
    if count_type === "norm"
        ct_obj = subset_count(sc_obj.normCount; genes = genes)
    elseif count_type === "raw"
        ct_obj = subset_count(sc_obj.rawCount; genes = genes)
    elseif count_type === "scale"
        ct_obj = subset_count(sc_obj.scaleCount; genes = genes)
    else
        println("count_type can only be \"raw\", \"norm\" or \"scale\"!")
    end
if isa(split_by, Nothing)
    all_df=DataFrame()
    for (i, gene) in enumerate(genes)
        gene_expr = subset_count(ct_obj; genes = [gene]).count_mtx
        gene_expr = vec(gene_expr)
        df = DataFrame()
        df.gene = gene_expr
        df.celltype=string.(sc_obj.metaData[!, cluster])
        avg_expr = DataFrames.combine(groupby(df, :celltype), :gene => mean => :avg_exp)
        perc_expr = DataFrames.combine(groupby(df, :celltype), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
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
        height= height, width=width
        )
else
    all_df=DataFrame()
    for (i, gene) in enumerate(genes)
        gene_expr = subset_count(ct_obj; genes = [gene]).count_mtx
        gene_expr = vec(gene_expr)
        df = DataFrame()
        df.gene = gene_expr
        df.celltype = string.(sc_obj.metaData[!, cluster])
        df.split_by = string.(sc_obj.metaData[!, split_by])
        avg_expr = DataFrames.combine(groupby(df, [:celltype, :split_by]), :gene => mean => :avg_exp)
        perc_expr = DataFrames.combine(groupby(df, [:celltype,:split_by]), :gene => function(x) countmap(x.>expr_cutoff)[:1]*100/length(x) end => :perc_exp)
        df_plt=innerjoin(avg_expr, perc_expr, on = [:celltype,:split_by])
        df_plt.gene .= gene
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
        height= height, width=width
        )        
    end
    return p
end

function violin_plot(sc_obj::Union{scRNAObject, VisiumObject, ImagingSpatialObject, CartanaObject, XeniumObject, scATACObject, MerfishObject, SlideseqObject, seqFishObject, STARmapObject,StereoSeqObject}, genes; 
    count_type::String ="norm", group_by::String = "cluster", 
    pt_size::Real =0.5, line_width::Real = 0, alpha::Real=1,
    height::Real = 800, width::Real = 500, do_legend::Bool = false,
    col_use::Union{Vector, Symbol, Nothing}=nothing)
    if count_type === "norm"
        ct_obj = subset_count(sc_obj.normCount; genes = genes)
    elseif count_type === "raw"
        ct_obj = subset_count(sc_obj.rawCount; genes = genes)
    elseif count_type === "scale"
        ct_obj = subset_count(sc_obj.scaleCount; genes = genes)
    else
        println("count_type can only be \"raw\", \"norm\" or \"scale\"!")
    end
    count_mat = ct_obj.count_mtx
    count_mat = count_mat'
    noise = randn(size(count_mat)[1]) ./ 1000
    count_mat = count_mat .+ noise
    count_mat = count_mat'
    count_mat = DataFrame(count_mat, :auto)
    count_mat.gene = ct_obj.gene_name
    count_mat = permutedims(count_mat, :gene)
    count_mat.cells = sc_obj.rawCount.cell_name
    count_mat.cluster = sc_obj.metaData[!, group_by]
    if isa(col_use, Nothing)
        col_use=:auto
    end
    gr(size=(width,height))
    l = @layout [a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z;aa;bb;cc;dd;ee;ff;gg;hh;ii;jj;kk;ll;mm;nn;oo;pp;qq;rr;ss;tt;uu;vv;ww;xx;yy;zz]
    l = l[1:length(genes)]
    p = []
    for (i, gene) in enumerate(genes)
        if i < length(genes)
        p1 = @df count_mat StatsPlots.violin(string.(:cluster), cols(Symbol.(gene)), ylabel=gene,
                group=string.(:cluster),linewidth=line_width, alpha=alpha, legend=do_legend, xaxis=nothing,grid = false,
                color_palette = col_use)
        p1 = @df count_mat dotplot!(string.(:cluster), cols(Symbol.(gene)), marker=(:black, stroke(1)), 
                    markersize =pt_size, legend=do_legend, xaxis=nothing,grid = false )
        else 
        p1 = @df count_mat StatsPlots.violin(string.(:cluster), cols(Symbol.(gene)), ylabel=gene,
                group=string.(:cluster),linewidth=line_width, alpha=alpha, legend=do_legend,grid = false, 
                color_palette = col_use)
        p1 = @df count_mat dotplot!(string.(:cluster), cols(Symbol.(gene)), marker=(:black, stroke(1)), 
                    markersize =pt_size, legend=do_legend,grid = false )
        end
        p = push!(p, p1)
    end
    return plot(p..., layout=l)
end