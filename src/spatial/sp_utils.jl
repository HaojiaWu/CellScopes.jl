function make_cell_proportion_df(cell_stat::DataFrame; 
        nfeatures::Int64=10, column::Union{Symbol, String}=:celltype, 
        time_point::Union{String, Nothing}=nothing)
    cell_stat=subset(cell_stat, :n_transcripts => ByRow(>=(nfeatures)))
    cell_counts=countmap(cell_stat[!, column])
    cell_counts=DataFrame(celltype=collect(keys(cell_counts)), count=collect(values(cell_counts)))
    total_cells=sum(cell_counts.count)
    cell_counts[!, :fraction].= 100*cell_counts.count/total_cells
    if time_point !== nothing
        cell_counts[!, :time] .= time_point
    end
    return cell_counts
end

function map_values(df::DataFrame,
        old_col::Union{Symbol, String},
        new_col::Union{Symbol, String},
        from::Union{Vector, String, Float64, Int64},
        to::Union{Vector, String, Float64, Int64})
    map_value=Dict(from .=> to)
    if isa(old_col, String)
        old_col = Symbol(old_col)
    end
    if isa(new_col, String)
        new_col = Symbol(new_col)
    end
    df = DataFrames.transform(df, old_col => ByRow(x -> map_value[x]) => new_col)
    return df
end

function reorder(df::DataFrame,
        col::Union{Symbol, String},
        order_by::Union{Vector, String})
    if isa(col, String)
        col = Symbol(col)
    end
    old_order=unique(df[!,col])
    order_index=collect(1:length(old_order))
    df = map_values(df, col, :order_index, order_by, order_index)
    sort!(df, :order_index)
    return df
end

function compute_pearson_cor(sp::Union{ImagingSpatialObject, CartanaObject, XeniumObject, MerfishObject, seqFishObject, STARmapObject, StereoSeqObject}, cluster1::Union{Symbol, String}, cluster2::Union{Symbol, String}; color_scheme::String="lightgreyred",reverse_color::Bool=false)
    df=sp.spmetaData.cell
    celltypes1=unique(df[!,cluster1])
    celltypes2=unique(df[!,cluster2])
    norm_cells=sp.normCounts
    df2=DataFrame([[] [] []],:auto)
    for i in celltypes1
        cols=DataFrame([[] [] []],:auto)
        for j in celltypes2
            x=filter(cluster1 => ==(i),df)
            cell1 = subset_count(norm_cells; cell_name=x.cell)
            cell1 = cell1.count_mtx
            dg = deepcopy(cell1)
            dg.rowmean .= mean(Array(dg), dims=2)
            x=dg.rowmean
            y=filter(cluster2 => ==(j),df)
            cell1 = subset_count(norm_cells; cell_name=y.cell)
            cell1 = cell1.count_mtx
            dg = deepcopy(cell1)
            dg.rowmean .= mean(Array(dg), dims=2)
            y=dg.rowmean
            cor_val=Statistics.cor(x,y)
            rows=DataFrame([string(i) string(j) cor_val],:auto)
            cols=[cols; rows]
        end
        df2=[df2;cols]
    end
    DataFrames.rename!(df2, [string(cluster1),string(cluster2),"cor"])
    df2.cor=float.(df2.cor)
    p=df2 |> @vlplot(:rect,  x=cluster1, y=cluster2, color={"cor:q",scale={scheme=color_scheme,reverse=reverse_color}})
    return p
end

function scan_cells(x, y, center_x, center_y, radius)
    sqrt((x - center_x)^2 + (y - center_y)^2)< radius
end

function compare_cell_distances(sp::Union{CartanaObject, XeniumObject,MerfishObject, seqFishObject, STARmapObject}, col::Union{String, Symbol}, target_cell::String, 
    cell1::String, cell2::String, radius::Union{Int64, Float64})
    coord_cells=sp.spmetaData.cell
    if isa(col, String)
        col=Symbol(col)
    end
    target=filter(col => ==(target_cell), coord_cells)
    target_cells=target.cell
    df = Array{Int64}(undef, 0, 2)
    n= length(target_cells)
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    Threads.@threads for i in target_cells
        ref_coord = filter(:cell => ==(i), coord_cells)
        center_x = ref_coord.x[1]
        center_y = ref_coord.y[1]
        within_range=filter([:x, :y]=>(a,b)->scan_cells(a, b, center_x, center_y, radius),coord_cells)
        cell_ct=countmap(within_range[!, col])
        if cell1 in keys(cell_ct)
            cell1_num=cell_ct[cell1]
        else
            cell1_num=0
        end
        if cell2 in keys(cell_ct)
            cell2_num=cell_ct[cell2]
        else
            cell2_num=0
        end
        row1=[cell1_num cell2_num]
        df = [df; row1]
        next!(p)
    end
    df=DataFrame(df, :auto)
    df = rename(df, :x1 => "cell1",:x2=>"cell2")
    total_cells=countmap(coord_cells[!, col])
    cell1_total=total_cells[cell1]
    cell2_total=total_cells[cell2]
    df.radius.=radius
    df.cell1_total.=cell1_total
    df.cell2_total.=cell2_total
    weighted = cell1_total/cell2_total
    df.cell1_density.=df.cell1/(radius^2*pi)
    df.cell2_density.=(weighted*df.cell2)/(radius^2*pi)
    return df
end

function coord_transform(point, degree;)
    degree = degree * pi/180
    x_rot = point[1]*cos(degree) + point[2]sin(degree)
    y_rot = -point[1]*sin(degree) + point[2]cos(degree)
    return (x_rot, y_rot)
end

function rotate_axis(sp::AbstractSpaObj, degree)
    coord_cells=sp.spmetaData.cell
    coord_molecules=sp.spmetaData.molecule
    new_x=[]
    new_y=[]
    for point in zip(coord_cells.x, coord_cells.y)
        new_point=coord_transform(point, degree)
        new_x = append!(new_x, [new_point[1]])
        new_y = append!(new_y, [new_point[2]])
    end
    coord_cells.x=float.(new_x)
    coord_cells.y=float.(new_y)
    new_x=[]
    new_y=[]
    for point in zip(coord_molecules.x, coord_molecules.y)
        new_point=coord_transform(point, degree)
        new_x = append!(new_x, [new_point[1]])
        new_y = append!(new_y, [new_point[2]])
    end
    coord_molecules.x=float.(new_x)
    coord_molecules.y=float.(new_y)
    sp.spmetaData.cell=coord_cells
    sp.spmetaData.molecule=coord_molecules
    return sp
end

function rotate_axis(df::DataFrame, x, y, degree)
    x_coord = df[!, x]
    y_coord = df[!, y]
    new_x=[]
    new_y=[]
    for point in zip(x_coord, y_coord)
        new_point=coord_transform(point, degree)
        new_x = append!(new_x, [new_point[1]])
        new_y = append!(new_y, [new_point[2]])
    end
    df.new_x=float.(new_x)
    df.new_y=float.(new_y)
    return df
end

function scale_data(X::Union{Vector{Int64}, Vector{Float64}}; 
    method::String="scanpy_scale", perc::Float64=0.02)
    if method ==="r-like-scale"
        scale_x = mapslices(Statistics.normalize!, X .- mean(X,dims=1), dims=1) 
        scale_x = scale_x * sqrt(size(X,1)-1)
    elseif method==="scanpy_scale"
        scale_x = clamp.(X, Statistics.quantile(X,perc), Statistics.quantile(X,1-perc))
    else
     error("method must be \"r-like-scale\" or \"scanpy_scale\"")
    end
    return scale_x
end

function unit_range_scale(data_input::Union{Vector{Int64}, Vector{Float64}})
    min_x=minimum(data_input)
    max_x=maximum(data_input)
    diff_x = max_x - min_x
    new_data = data_input .- min_x
    new_data = new_data ./ diff_x
    return new_data
end

function split_field(df::DataFrame, n_fields_x::Int64, n_fields_y::Int64)
    x_min = minimum(df.x)-1
    x_max = maximum(df.x)+1
    y_min = minimum(df.y)-1
    y_max = maximum(df.y)+1
    x_seg = (x_max-x_min)/n_fields_x
    y_seg = (y_max-y_min)/n_fields_y
    pts2=[]
    centroids2=()
    for i in 1:n_fields_x
      pts=[]
      centroids=()
        for j in 1:n_fields_y
          p1=MK.Point2f0(x_min+x_seg * (i-1), y_min+y_seg * (j-1))
          p2=MK.Point2f0(x_min+x_seg * i, y_min+y_seg * (j-1))
          p3=MK.Point2f0(x_min+x_seg * i, y_min+y_seg * j)
          p4=MK.Point2f0(x_min+x_seg * (i-1), y_min+y_seg * j)
          centroid=(x_min + x_seg * i - x_seg/2, y_min + y_seg * j - y_seg/2)
          rect_pts=[p1,p2,p3,p4]
          pts=append!(pts, [rect_pts])
          centroids=[centroids...,centroid]
        end
      pts2=[pts2 ; pts]
      centroids2=[centroids2;centroids]
    end
  deleteat!(centroids2,1)
  return pts2, centroids2
end

function slope2deg(slope::Float64)
    if slope >=0 
        degree=rad2deg(atan(slope))
    else
        degree=180-abs(rad2deg(atan(slope)))
    end
    return degree
end

function subset_fov(sp::Union{ImagingSpatialObject, CartanaObject, XeniumObject,MerfishObject, seqFishObject, STARmapObject, StereoSeqObject, VisiumObject}, fov::Vector{Int64}, n_fields_x::Int64, n_fields_y::Int64)
    if isa(sp, VisiumObject)
        df = deepcopy(sp.spmetaData)
        rename!(df, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, :x, :y])
    else
        df = deepcopy(sp.spmetaData.cell)
    end
    pts, centroids = split_field(df, n_fields_x, n_fields_y)
    pts_sub=pts[fov]
    min_pt=minimum(minimum(pts_sub))
    max_pt=maximum(maximum(pts_sub))
    min_pt=convert(Vector{Float32}, min_pt)
    max_pt=convert(Vector{Float32}, max_pt)
    df_sub=filter([:x, :y] => (x, y) -> x > min_pt[1] && x < max_pt[1] && y > min_pt[2] && y < max_pt[2], df)
    return df_sub      
end

function get_slope(pt1, pt2)
    slope= (pt1[2]-pt2[2])/(pt1[1]-pt2[1])
    return slope
end

function compute_new_coord(df, pt, center; span=150)
    slope=get_slope(pt, center)
    degree=slope2deg(slope)
    pt_new0=[pt[1]-center[1], pt[2]-center[2]]
    pt_new=coord_transform(pt_new0, degree)
    x_new =[]
    y_new=[]
    for (x, y) in zip(df.x, df.y)
        pt1=(x,y)
        pt2=(pt1[1]-center[1], pt1[2]-center[2])
        pt3 = coord_transform(pt2, degree)
        x_new = append!(x_new, pt3[1])
        y_new = append!(y_new, pt3[2])
    end
    new_coord=DataFrame(x_new=x_new, y_new=y_new)
    new_coord2 = new_coord[findall((-span) .< new_coord.y_new .< span),:]
    if pt_new[1]<0
       length_x=minimum(new_coord2.x_new)
    else 
       length_x=maximum(new_coord2.x_new)
    end
    depth=abs((pt_new[1]-length_x)/length_x)
    angle=atan(pt_new0[2],pt_new0[1])
    return depth, angle
end

function compute_kidney_coordinates(sp::Union{ImagingSpatialObject, CartanaObject, XeniumObject, MerfishObject, seqFishObject, STARmapObject, StereoSeqObject}, center)
    df = sp.spmetaData.cell
    kid_depth=[]
    kid_angle=[]
    n=length(df.x)
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    Threads.@threads for i in 1:n
        pt=(df.x[i], df.y[i])
        depth, angle = compute_new_coord(df, pt, center)
        kid_depth=append!(kid_depth, depth)
        kid_angle=append!(kid_angle, angle)
        next!(p)
    end
    df[!,:depth]=kid_depth
    df[!,:angle]=kid_angle
    sp.spmetaData.cell=df
    molecules=sp.spmetaData.molecule
    cells2=df.cell
    molecules=filter(:cell=> ∈(Set(cells2)), molecules)
    from=df.cell
    to=df.depth
    molecules_df=map_values(molecules, :cell, :depth,from, to)
    sp.spmetaData.molecule=molecules_df
    return sp
end

function compute_kidney_coordinates(df::DataFrame, x_col, y_col, center)
    kid_depth=[]
    kid_angle=[]
    kid_length=[]
    n=length(df[!, x_col])
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    Threads.@threads for i in 1:n
        pt=(df[!, x_col][i], df[!, y_col][i])
        depth, angle, total_length = compute_new_coord(df,x_col ,y_col,pt, center; span1=100)
        kid_depth=append!(kid_depth, depth)
        kid_angle=append!(kid_angle, angle)
        kid_length=append!(kid_length, total_length)
        next!(p)
    end
    kid_length=float.(kid_length)
    kid_angle=float.(kid_angle)
    kid_depth=float.(kid_depth)
    df[!,:depth]=kid_depth
    df[!,:angle]=kid_angle
    df[!,:total_length]=kid_length
    return df
end

function point_dist(pt1::Union{Tuple, Vector}, pt2::Union{Tuple, Vector})
    pt_dist = sqrt((pt1[1]-pt2[1])^2 + (pt1[2]-pt2[2])^2)
    return pt_dist
end

function align_coordinates(df::DataFrame, x_col::Union{String, Symbol}, y_col::Union{String, Symbol})
    x_coord = df[!,x_col]
    y_coord = df[!,y_col]
    center_x = mean(x_coord)
    center_y = mean(y_coord)
    all_points=[]
    for point in zip(x_coord, y_coord)
        all_points=append!(all_points, [point])
    end
    dist0=point_dist(all_points[1], all_points[2])
    pos_max=[]
    n = length(all_points)-1
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    Threads.@threads for i in 1:(length(all_points)-1)
        pt1 = all_points[i]
        for j in (i+1):length(all_points)
            pt2 = all_points[j]
            pt_dist = point_dist(pt1, pt2)
            if pt_dist > dist0
                dist0 = pt_dist
                pos_max = [i, j]
            else
                continue
            end
        end
    next!(p)
    end
    pt1_max = all_points[pos_max[1]]
    pt2_max = all_points[pos_max[2]]
    pt_x = [pt1_max[1], center_x, pt2_max[1]]
    pt_y = [pt1_max[2], center_y, pt2_max[2]]
    quadfit = Polynomials.fit(pt_x,pt_y,1)
    slope = quadfit[1]
    degree = slope2deg(slope)
    df_new = rotate_axis(df, x_col, y_col, degree)
    return df_new
end

function read_seurat_data(file_name; file_format="rda")
    if file_format === "rda"
        rbase = rimport("base")
        rbase.rm(list =rbase.ls())
        rbase.load(file_name)
        seu = rbase.get(rbase.ls())
        return seu
    elseif file_format === "rds"
        seu = R"readRDS"(file_name)
        return seu
    else
        error("file_format can only be \"rda\" or \"rds\"")
    end
end

function read_scanpy_data(file_name)
    sc = pyimport("scanpy")
    adata = sc.read(file_name)
    return adata
end

function visium_actual_dist(spot_num)
    vs_dist = 100 * spot_num + 55    #Spot diameter + c-c distance
    return vs_dist
end

function visium_unit_radius(spot_num)
    vs_dist = visium_actual_dist(spot_num)
    r_unit = 55/2/vs_dist
    return r_unit
end

function visium_deconvolution(vs::Union{VisiumObject, SlideseqObject}, sp::Union{ImagingSpatialObject, CartanaObject, XeniumObject, MerfishObject, seqFishObject, STARmapObject, StereoSeqObject}, spot_r::Union{Int64, Float64};
    vscell_col = "cell" , spcluster_col="celltype", vs_x = "new_x", vs_y = "new_y", 
    sp_x = "new_x", sp_y = "new_y")
    vs_cells = deepcopy(vs.cells)
    sp_cells = deepcopy(sp.spmetaData.cell)
    target_cells = vs_cells.cell
    celltypes = string.(keys(countmap(sp_cells[!, spcluster_col])))
    df1 = DataFrame()
    for i in target_cells
        ref_coord = filter(vscell_col => ==(i), vs_cells)
        center_x = ref_coord[!, vs_x][1]
        center_y = ref_coord[!, vs_y][1]
        within_range=filter([sp_x, sp_y] => (a,b) -> scan_cells(a, b, center_x, center_y, spot_r), sp_cells)
        cell_ct = countmap(within_range[!, spcluster_col])
        total_num = sum(values(cell_ct))
        if cell_ct == Dict{String, Int64}()
            cell_ct = Dict(celltypes .=> 0)
        end
        df2 = DataFrame(cell_ct...)
        for j in celltypes
            if !(j in keys(cell_ct))
                cell1_num=0
                df2[!, j] .= cell1_num
            end
        end
        df2=df2[!, celltypes]
        df2.cell .= i
        df2.total .= total_num
        df1 = [df1; df2]
    end
    new_df =innerjoin(vs_cells, df1, on = :cell)
    return new_df
end

function split_spatial(n_bin::Int64)
    n_seg = 1.0/n_bin
    all_segs=[]
    for i in 1:n_bin
        seg_low = n_seg * (i-1)
        seg_high = n_seg * i
        seg_i =[seg_low,seg_high]
        all_segs=append!(all_segs, [seg_i])
    end
    return all_segs
end

function bin_gene_spatial(sp::Union{ImagingSpatialObject, CartanaObject, XeniumObject,MerfishObject, seqFishObject, STARmapObject}, n_bin::Int64; 
    celltype::Union{String, Int64, Nothing}=nothing,
    assay_use = "measured", imp_type = "SpaGE", genes = nothing)
    cells=deepcopy(sp.spmetaData.cell)
    if celltype !== nothing
        cells = filter(:celltype => ==(celltype), cells)
    end
    all_segs = split_spatial(n_bin)
    n_seg = 1.0/n_bin
    new_df = DataFrame()
    for i in 1:length(all_segs)
        seg = all_segs[i]
        cell1 = cells[findall(seg[1] .< cells.depth .< seg[2]),:]
        cell1.bin .= n_seg * i
        new_df = [new_df; cell1]
    end
    if assay_use === "measured"
        gene_expr = deepcopy(sp.normCount)
    elseif assay_use === "predicted"
        if imp_type === "tangram"
            gene_expr = sp.imputeData.tgCount
        elseif imp_type === "SpaGE"
            gene_expr = sp.imputeData.spageCount
        elseif imp_type === "gimVI"
            gene_expr = sp.imputeData.gimviCount
        else
            error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
        end
    else
        error("assay_use can only be \"measured\" or \"predicted\"")
    end
    if isa(genes, Nothing)
        all_genes = gene_expr.gene_name
    else
        all_genes = genes
    end
    gene_expr = subset_count(gene_expr; genes = all_genes)
    gene_expr = ctobj_to_df(gene_expr)
    df_proc = innerjoin(gene_expr, new_df, on = :cell)
    new_df2 = DataFrame()
    for gene in all_genes
        avg_expr=DataFrames.combine(groupby(df_proc, :bin), gene => mean => :avg_exp)
        avg_expr.gene .= gene
        new_df2 = [new_df2; avg_expr]
    end
    new_df2.gene = string.(new_df2.gene);
    return new_df2
end

#= This is deprecated.
function pd_to_df(df_pd)
    colnames = map(Symbol, df_pd.columns)
    df = DataFrame(Any[Array(df_pd[c].values) for c in colnames], colnames)
end
=#

function pd_to_df(pydf)
    cols = pydf.columns.to_list()
    jldf = DataFrame(Any[collect(values(pydf[c])) for c in cols], map(Symbol, cols))
    return jldf
end

function ctobj_to_df(ct_obj::AbstractCount)
    count_mat = ct_obj.count_mtx'
    count_df = DataFrame(count_mat, :auto)
    rename!(count_df, ct_obj.gene_name)
    count_df = hcat(ct_obj.cell_name, count_df)
    rename!(count_df, :x1 => :cell)
    return count_df
end

function gunzip(fname) ## function from https://github.com/JuliaSparse/MatrixMarket.jl/blob/570c6d4e443d464605264b30c76b32bec102b3b8/test/dl-matrixmarket.jl
    destname, ext = splitext(fname)
    if ext != ".gz"
        error("gunzip: $fname: unknown suffix -- ignored")
    end
    open(destname, "w") do f
        GZip.open(fname) do g
            write(f, read(g, String))
        end
    end
    destname
end

function log_norm_spatial(x, median_value)
    return log.(((x ./ sum(x)) .* median_value) .+ 1)
end

function log_norm_cpm(x)
    return log.(((x ./ sum(x)) .* 1000000) .+ 1)
end

function from_scanpy(adata::Union{String, PyObject}; data_type = "scRNA", 
    tech::Union{String, Nothing}=nothing,
    sp_coord_name="spatial", anno="leiden")
    sc = pyimport("scanpy")
    if isa(adata, String)
        adata = sc.read(adata)
    end
    genes = convert(Vector{String},adata.var_names)
    cells = convert(Vector{String}, adata.obs_names)
    umap = get(adata.obsm,"X_umap")
    meta = pd_to_df(adata.obs)
    meta.cell = cells
    if data_type == "scRNA"
        scale_count = convert(SparseMatrixCSC,adata.X)
        scale_count = scale_count'
        sim_count = sprand(size(scale_count)[1],size(scale_count)[2],0.1) ## simulate a count matrix but it is not going to use
        raw_count = RawCountObject(sim_count, cells, genes)     
        cs_obj = scRNAObject(raw_count;meta_data=meta)
    elseif data_type == "spatial"
        if isa(tech, Nothing)
            error("Please specify the spatial technology names in the 'tech' parameter. It can be 'xenium' or 'visium'.")
        elseif tech == "xenium"
            scale_count = convert(SparseMatrixCSC, adata.X.transpose())
            sim_count = sprand(size(scale_count)[1],size(scale_count)[2],0.1) ## simulate a count matrix but it is not going to use        
            raw_count = RawCountObject(sim_count, cells, genes)
            sp_coord = adata.obsm[sp_coord_name]
            sp_coord = DataFrame(sp_coord, :auto)
            rename!(sp_coord, [:x, :y])
            sp_coord.cell = convert(Vector{String}, adata.obs_names)
            sp_coord[!,anno] = meta[!, anno]
            cs_obj = ImagingSpatialObject(sp_coord, sp_coord,raw_count; meta_data=meta)
        elseif tech == "visium"
            scale_count = convert(SparseMatrixCSC, adata.X.transpose())
            sim_count = sprand(size(scale_count)[1],size(scale_count)[2],0.1) ## simulate a count matrix but it is not going to use
            raw_count = RawCountObject(sim_count, cells, genes)
            sp_coord = adata.obsm[sp_coord_name]
            sp_coord = DataFrame(sp_coord, :auto)
            rename!(sp_coord, [:x, :y])
            sp_coord.cell = cells
            sp_coord[!,anno] = meta[!, anno]
            rename!(sp_coord, [:cell, :y, :x] .=> [:barcode, :pxl_col_in_fullres, :pxl_row_in_fullres])
            img_high = adata.uns["spatial"][collect(keys(adata.uns["spatial"]))[1]]["images"]["hires"]
            img_high = [RGB(img_high[i, j, 1], img_high[i, j, 2], img_high[i, j, 3]) for i in 1:size(img_high, 1), j in 1:size(img_high, 2)]
            img_high = img_high'
            img_high = convert(Matrix{RGB{N0f8}}, img_high)
            img_low = adata.uns["spatial"][collect(keys(adata.uns["spatial"]))[1]]["images"]["lowres"]
            img_low = [RGB(img_low[i, j, 1], img_low[i, j, 2], img_low[i, j, 3]) for i in 1:size(img_low, 1), j in 1:size(img_low, 2)]
            img_low = img_low'
            img_low = convert(Matrix{RGB{N0f8}}, img_low)
            json_data = adata.uns["spatial"][collect(keys(adata.uns["spatial"]))[1]]["scalefactors"]
            json_data = convert(Dict{String, Any}, json_data)
            cs_obj = VisiumObject(raw_count)
            image_obj = VisiumImgObject(img_high, img_low, nothing, nothing, nothing, json_data)
            cs_obj.metaData = meta
            cs_obj.imageData = image_obj
            cs_obj.spmetaData = sp_coord
        else
            error("tech can only be visum or xenium for now!")
        end
    else
        error("Currently only support scRNA and spatial")
    end
    cs_obj = normalize_object(cs_obj)
    cs_obj = scale_object(cs_obj)
    cs_obj.rawCount.count_mtx = scale_count
    cs_obj.normCount.count_mtx = scale_count
    cs_obj.scaleCount.count_mtx = scale_count ### replace all count with the scale_count from scanpy
    umap = UMAPObject(umap, 
        "UMAP", 
        size(umap)[2],
        nothing,
        nothing, 
        nothing,
        nothing,
        nothing
        )
    cs_obj.dimReduction = ReductionObject(nothing, nothing, umap)
    clusters = DataFrame(:cell=>cells, :cluster => meta[!, anno])
    cs_obj.clustData = ClusteringObject(clusters, 
        nothing, 
        nothing, nothing, nothing )
    return cs_obj
end

function from_seurat(seurat_file; data_type::String = "scRNA", 
                tech::Union{String, Nothing} = nothing,  
                anno::String = "cluster", 
                assay::String = "RNA", 
                version::String = "v5")
    rbase = rimport("base")
    seu = rimport("Seurat")
    seu_obj = rbase.readRDS(seurat_file)
    meta = seu_obj["meta.data"]
    meta = rcopy(meta)
    clusters = rcopy(seu_obj["active.ident"])
    meta[!, anno] = clusters
    if version=="v5"
        counts = rcopy(seu_obj["assays"][assay]["layers"]["counts"])
    else
        counts = rcopy(seu_obj["assays"][assay]["counts"])
    end
    counts = sparse_r_to_jl(counts)
    umap = rcopy(seu_obj["reductions"]["umap"]["cell.embeddings"])
    cells = rcopy(rbase.colnames(seu_obj))
    genes = rcopy(rbase.rownames(seu_obj))
    umap = UMAPObject(umap, 
        "UMAP", 
        size(umap)[2],
        nothing,
        nothing, 
        nothing,
        nothing,
        nothing
    )
    raw_count = RawCountObject(counts, cells, genes)
    if data_type == "scRNA"
        cs_obj = scRNAObject(raw_count; meta_data=meta)
        if version=="v5"
            norm_count = rcopy(seu_obj["assays"][assay]["layers"]["data"])
        else
            norm_count = rcopy(seu_obj["assays"][assay]["data"])
        end
        norm_count = sparse_r_to_jl(norm_count)
        cs_obj = normalize_object(cs_obj)
        cs_obj.normCount.count_mtx = norm_count
    elseif data_type == "spatial"
        if isa(tech, Nothing)
            error("Please specify the spatial technology names in the 'tech' parameter. It can be 'xenium' or 'visium'.")
        elseif tech == "xenium"
            x_coords = Float64[]
            y_coords = Float64[]
            gene_symbols = String[]
            for gene in genes
                gene_coord = rcopy(seu_obj["images"]["fov"]["molecules"]["molecules"][gene]["coords"])
                append!(x_coords, gene_coord[:, 1])
                append!(y_coords, gene_coord[:, 2])
                append!(gene_symbols, fill(gene, size(gene_coord, 1)))
            end
            molecules = DataFrame(x = x_coords, y = y_coords, gene = gene_symbols)
            x_coords = Float64[]
            y_coords = Float64[]
            cell_symbols = String[]
            for i in 1:length(cells)
                poly = rcopy(seu_obj["images"]["fov"]["boundaries"]["segmentation"]["polygons"][i]["Polygons"][1]["coords"])
                append!(x_coords, poly[:, 1])
                append!(y_coords, poly[:, 2])
                append!(cell_symbols, fill(cells[i], size(poly, 1)))
            end
            seg = DataFrame(x = x_coords, y = y_coords, cell_id = cell_symbols)
            grouped = groupby(seg, :cell_id)
            cell_ids = unique(seg.cell_id)
            poly = Vector{Matrix{Float64}}(undef, length(cell_ids))
            n = length(cell_ids)
            println("Formatting cell polygons...")
            p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
            for idx in 1:length(cell_ids)
                cell_data = grouped[idx]
                cell_1 = Matrix(cell_data[!, 1:2])
                poly[idx] = cell_1
                next!(p)
            end
            cell_coords = rcopy(seu.GetTissueCoordinates(seu_obj))
            meta.cell = cell_coords.cell
            cell_coords[!, anno] = rcopy(seu_obj["active.ident"])
            cs_obj = XeniumObject(molecules, cell_coords, raw_count, poly, umap; meta_data=meta)
        elseif tech == "visium"
            cs_obj = VisiumObject(raw_count)
            positions = rcopy(seu.GetTissueCoordinates(seu_obj))
            positions.barcode = cells
            positions.cluster = clusters
            rename!(positions, :imagecol => :pxl_row_in_fullres, :imagerow =>:pxl_col_in_fullres)
            cs_obj.spmetaData = positions
            keys_scale = ["tissue_lowres_scalef", "tissue_hires_scalef"]
            values_scale = [1 , 1]
            json_data = Dict{String, Any}(keys_scale .=> values_scale)
            img = rcopy(seu_obj["images"][1]["image"])
            img = [RGB(img[i, j, 1], img[i, j, 2], img[i, j, 3]) for i in 1:size(img, 1), j in 1:size(img, 2)]
            img = img'
            img = convert(Matrix{RGB{N0f8}}, img)
            image_obj = VisiumImgObject(img, img, nothing, nothing, nothing, json_data)
            cs_obj.metaData = meta
            cs_obj.imageData = image_obj
            cs_obj = normalize_object(cs_obj)
        else
            error("tech can only be visum or xenium for now!")
        end
    else
        error("data_type can only be scRNA or spatial")
    end
    cs_obj.dimReduction = ReductionObject(nothing, nothing, umap)
    clusters = DataFrame(:cell=>cells, :cluster => clusters)
    cs_obj.clustData = ClusteringObject(clusters, 
        nothing, 
        nothing, nothing, nothing )
    return cs_obj
end

function spatial_range(sp; x_col::String="x", y_col::String="y")
    if isa(sp, VisiumObject)
        coord_cell = deepcopy(sp.spmetaData)
        x_col = Symbol(x_col)
        y_col = Symbol(y_col)
        rename!(coord_cell, [:barcode, :pxl_row_in_fullres, :pxl_col_in_fullres] .=> [:cell, x_col, y_col])
    else
        coord_cell = deepcopy(sp.spmetaData.cell)
    end
        x_min = minimum(coord_cell[!,x_col])
        x_max = maximum(coord_cell[!,x_col])
        y_min = minimum(coord_cell[!,y_col])
        y_max = maximum(coord_cell[!,y_col])
    return [(x_min, x_max); (y_min, y_max)]
end

get_xn_sf = function(sp::XeniumObject;adjust_coord_to_img="auto", x_col="x", y_col="y")
    scale_values = Dict(
                        "level1" => 0.2125 / 0.2125,
                        "level2" => 0.4250 / 0.2125,
                        "level3" => 0.8500 / 0.2125,
                        "level4" => 1.7000 / 0.2125,
                        "level5" => 3.4000 / 0.2125,
                        "level6" => 6.8000 / 0.2125,
                        "level7" => 13.6000 / 0.2125,
                        "level8" => 27.2000 / 0.2125
                        )
    scale_value = get(scale_values, adjust_coord_to_img, 0)
    scale_x = scale_y = scale_value
    if scale_value ==0
        img = deepcopy(sp.imageData)
        df_plt = deepcopy(sp.spmetaData.cell)
        scale_x = maximum(df_plt[!, x_col]) / size(img)[1]
        scale_y = maximum(df_plt[!, y_col]) / size(img)[2]
    end
    return (scale_x, scale_y)
end

get_vs_sf = function(sp::VisiumObject;img_res="high")
        if img_res == "high"
            scale_factor = sp.imageData.jsonParameters["tissue_hires_scalef"]
        elseif img_res == "low"
            scale_factor = sp.imageData.jsonParameters["tissue_lowres_scalef"]
        elseif img_res == "full"
            dim_full = size(sp.imageData.fullresImage)
            dim_high = size(sp.imageData.highresImage)
            x_ratio = dim_full[1]/dim_high[1]
            y_ratio = dim_full[2]/dim_high[2]
            scale_factor = sp.imageData.jsonParameters["tissue_hires_scalef"]
            scale_factor = scale_factor * (x_ratio + y_ratio)/2    
        else
            error("img_res can only be \"high\", \"low\" or \"full\"!")
        end
        return scale_factor
end

function generate_alternating_pattern(n, m; is_reverse::Bool = true)
    ascending = 0:(n-1)
    descending = (n-1):-1:0
    result = Int[]
    for i = 1:m
        if is_reverse
            if isodd(i)
                append!(result, descending)
            else
                append!(result, ascending)
            end
        else
            if isodd(i)
                append!(result, ascending)
            else
                append!(result, descending)
            end
        end
    end
    return result
end

function generate_repeating_pattern(n, m)
    series_num = (m-1):-1:0
    result = repeat(series_num, inner=n)
    return result
end

function invert_y_axis(y_values::Union{Vector{Int64}, Vector{Float64}})
    min_y, max_y = minimum(y_values), maximum(y_values)
    inverted_y = max_y + min_y .- y_values
    return inverted_y
end

function make_ct_from_tx(tx_df; cell_col="CellId", gene_col="target")
    if isa(cell_col, String)
        cell_col = Symbol(cell_col)
    end
    if isa(gene_col, String)
        gene_col = Symbol(gene_col)
    end
    tx_df[!, :count] = ones(Int, nrow(tx_df))
    grouped_data = groupby(tx_df, [cell_col, gene_col])
    combined_data = combine(grouped_data, :count => sum)
    gene_count = unstack(combined_data, gene_col, cell_col, :count_sum, fill=0)
    rename!(gene_count, gene_col => :gene)
    return gene_count
end

function make_ct_from_df(df; gene_col = "gene")
    df[!, gene_col] = String.(df[!, gene_col])
    test_val = df[!, gene_col][1]
    if !isa(test_val, String)
        error("The first column of the dataframe must be set to gene column!")
    end
    gene_name = df[!, gene_col]
    cell_name = names(df)[2:end]
    count_df = df[!, 2:end]
    count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
    raw_count = RawCountObject(count_df, cell_name, gene_name)
    return raw_count
end


