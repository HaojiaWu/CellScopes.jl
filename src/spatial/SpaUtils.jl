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

function mapvalues(df::DataFrame,
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
    df = mapvalues(df, col, :order_index, order_by, order_index)
    sort!(df, :order_index)
    return df
end

function compute_pearson_cor(sp::CartanaObject, cluster1::Union{Symbol, String}, cluster2::Union{Symbol, String}; color_scheme::String="lightgreyred",reverse_color::Bool=false)
    df=sp.cells
    celltypes1=unique(df[!,cluster1])
    celltypes2=unique(df[!,cluster2])
    norm_cells=sp.norm_counts
    df2=DataFrame([[] [] []],:auto)
    for i in celltypes1
        cols=DataFrame([[] [] []],:auto)
        for j in celltypes2
            x=filter(cluster1 => x-> x == i,df)
            cell1=norm_cells[!, x.cell2]
            dg = deepcopy(cell1)
            dg.rowmean .= mean(Array(dg), dims=2)
            x=dg.rowmean;
            y=filter(cluster2 => y-> y == j,df)
            cell1=norm_cells[!, y.cell2]
            dg = deepcopy(cell1)
            dg.rowmean .= mean(Array(dg), dims=2)
            y=dg.rowmean;
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

function compare_cell_distances(sp::CartanaObject,col::Union{String, Symbol}, target_cell::String, 
    cell1::String, cell2::String, radius::Union{Int64, Float64})
    coord_cells=sp.cells
    if isa(col, String)
        col=Symbol(col)
    end
    target=filter(col => x -> x ==target_cell, coord_cells)
    target_cells=target.cell2;
    df = Array{Int64}(undef, 0, 2)
    n= length(target_cells)
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    Threads.@threads for i in target_cells
        ref_coord = filter(:cell2 => x-> x== i, coord_cells)
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
    coord_cells=sp.cells
    coord_molecules=sp.molecules
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
    sp.cells=coord_cells
    sp.molecules=coord_molecules
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

function subset_fov(sp::CartanaObject, fov::Vector{Int64}, n_fields_x::Int64, n_fields_y::Int64)
    df=sp.cells
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
    new_coord=DataFrame(x_new=x_new, y_new=y_new);
    new_coord2=filter(:y_new => y-> (-span) < y < span, new_coord)
    if pt_new[1]<0
       length_x=minimum(new_coord2.x_new)
    else 
       length_x=maximum(new_coord2.x_new)
    end
    depth=abs((pt_new[1]-length_x)/length_x)
    angle=atan(pt_new0[2],pt_new0[1])
    return depth, angle
end

function compute_kidney_coordinates(sp::CartanaObject, center)
    df = sp.cells
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
    sp.cells=df
    molecules=sp.molecules
    cells2=df.cell2
    molecules=filter(:cell2=> x-> x in cells2, molecules)
    from=df.cell2
    to=df.depth
    molecules_df=mapvalues(molecules, :cell2, :depth,from, to)
    sp.molecules=molecules_df
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

function visium_deconvolution(vs::VisiumObject,sp::CartanaObject, spot_r::Union{Int64, Float64};
    vscell_col = "cell" , spcluster_col="celltype", vs_x = "new_x", vs_y = "new_y", 
    sp_x = "new_x", sp_y = "new_y")
    vs_cells = deepcopy(vs.cells)
    sp_cells = deepcopy(sp.cells)
    target_cells = vs_cells.cell
    celltypes = string.(keys(countmap(sp_cells[!, spcluster_col])))
    df1 = DataFrame()
    for i in target_cells
        ref_coord = filter(vscell_col => x-> x == i, vs_cells)
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

function bin_gene_spatial(sp::CartanaObject, n_bin::Int64; celltype::Union{String, Int64, Nothing}=nothing)
    cells=deepcopy(sp.cells)
    if celltype !== nothing
        cells = filter(:celltype => x -> x == celltype, cells)
    end
    all_segs = split_spatial(n_bin)
    n_seg = 1.0/n_bin
    new_df = DataFrame()
    for i in 1:length(all_segs)
        seg = all_segs[i]
        cell1 = filter(:depth => x -> seg[1]< x <= seg[2], cells)
        cell1.bin .= n_seg * i
        new_df = [new_df; cell1]
    end
    gene_expr = deepcopy(sp.norm_counts)
    all_genes = gene_expr.gene
    gene_expr = permutedims(gene_expr, 1)
    rename!(gene_expr, :gene => :cell2)
    df_proc = innerjoin(gene_expr, new_df, on = :cell2)
    new_df2 = DataFrame()
    for gene in all_genes
        avg_expr=combine(groupby(df_proc, :bin), gene => mean => :avg_exp)
        avg_expr.gene .= gene
        new_df2 = [new_df2; avg_expr]
    end
    new_df2.gene = string.(new_df2.gene);
    return new_df2
end

function pd_to_df(df_pd)
    colnames = map(Symbol, df_pd.columns)
    df = DataFrame(Any[Array(df_pd[c].values) for c in colnames], colnames)
end

