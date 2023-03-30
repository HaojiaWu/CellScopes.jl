function polygons_cell_mapping(sp::Union{CartanaObject, XeniumObject}; anno::Union{String, Symbol}="cluster")
    polygons=sp.polygonData
    center_df=DataFrame()
    cell=Int[]
    X=Float64[]
    Y=Float64[]
    for (i,p) in enumerate(polygons)
        x=mean(p[:,1])
        y=mean(p[:,2])
        push!(cell,i)
        push!(X,x)
        push!(Y,y)
    end
    center_df.cell=cell
    center_df.X=X
    center_df.Y=Y
    cell_coord=sp.spmetaData.cell
    poly_mtx = Matrix(center_df[:,2:3])
    ref_mtx = Matrix(cell_coord[:,2:3])
    ref_coord=[]
    for i in eachrow(ref_mtx)
        ref_coord=push!(ref_coord,i)
    end
    query_coord=[]
    for i in eachrow(poly_mtx)
        query_coord=push!(query_coord,i)
    end
    n=length(query_coord)
    query_dist=[]
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for j in 1:length(query_coord)
        ref_dist=[Distances.euclidean(query_coord[j],ref_coord[i]) for i in 1:length(ref_coord)]
        query_dist=push!(query_dist,findmin(ref_dist)[2])
        next!(p)
    end
    println("Polygons have been mapped to the cells!")
    my_annot=DataFrame(polygon_number=collect(1:length(query_dist)), mapped_cell=query_dist);
    from=sp.spmetaData.cell.cell
    to=sp.spmetaData.cell[!, anno]
    prefix = split(sp.spmetaData.cell.cell[1], "_")
    if length(prefix) > 1
        my_annot.mapped_cell = prefix[1] * "_" .* string.(my_annot.mapped_cell);
    end
    anno=map_values(my_annot, :mapped_cell, :cluster,from, to)
    sp.spmetaData.polygon=anno
    return sp
end

function generate_polygon_counts(sp::AbstractSpaObj)
    coord_molecules = deepcopy(sp.spmetaData.molecule)
    if isa(sp.spmetaData.polygon, Nothing)
        error("Please run polygons_cell_mapping first!")
    end
    anno = deepcopy(sp.spmetaData.polygon)
    coord_molecules.cell = string.(coord_molecules.cell)
    anno = rename!(anno, :polygon_number => :number)
    anno = rename!(anno, :mapped_cell => :cell)
    anno.cell = string.(anno.cell)
    join_df = innerjoin(coord_molecules[!,[:gene,:cell]],anno[!,[:number,:cell]], on=:cell)
    join_df.gene = string.(join_df.gene)
    gdf = groupby(join_df, :gene)
    new_df = DataFrame()
    for (i, df) in enumerate(gdf)
        map_dict = countmap(df.number)
        anno1 = DataFrame(cell=collect(keys(map_dict)),count=collect(values(map_dict)))
        anno1.gene .= gdf[i].gene[1]
        new_df = [new_df;anno1]
    end
    final_df = unstack(new_df, :cell, :gene, :count)
    final_df .= ifelse.(isequal.(final_df, missing), 0, final_df)
    final_df = mapcols(ByRow(Int64), final_df);
    sort!(final_df, :cell)
    prefix = split(sp.spmetaData.cell.cell[1], "_")
    if length(prefix) > 1
        final_df.cell = prefix[1] * "_" .* string.(final_df.cell)
    end
    cellnames = final_df.cell
    genenames = names(final_df)[2:end]
    final_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(final_df[!,2:end])')
    raw_count = RawCountObject(final_df, cellnames, genenames)
    norm_count = normalize_object(raw_count)
    sp.polyCount = raw_count
    sp.polynormCount = norm_count
    println("Polygon counts was normalized!")
    return sp
end
