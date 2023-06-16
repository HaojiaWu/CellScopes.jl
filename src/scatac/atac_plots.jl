function gene_activity_plot(atac_obj::scATACObject, genes; dim_type::String = "umap", count_type = "norm",x_lims=nothing, y_lims=nothing, marker_size=4, order=true,
    color_keys::Union{Vector{String}, Tuple{String,String,String}}=("black","yellow","red"), do_dimname::Bool=false,
        split_by::Union{String, Symbol, Nothing}=nothing, titlesize::Int64 = 24, height::Real = 500, width::Real = 500)
        dim_data, x_col, y_col = get_dim_data(atac_obj.dimReduction, dim_type)
        ct_obj = subset_count(atac_obj.activityData; genes = genes)
        count_mat = ct_obj.count_mtx
        count_mat = DataFrame(count_mat, :auto)
        count_mat.gene = ct_obj.gene_name
        count_mat = permutedims(count_mat, :gene)
        count_mat.cells = atac_obj.rawCount.cell_name
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
                    n_col1 = 2*(i-3*(n_rows-1))-1
                    n_col2 = 2*(i-3*(n_rows-1))
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
            group_arr = string.(atac_obj.metaData[!, split_by])
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