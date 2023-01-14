function polygons_cell_mapping(sp::CartanaObject)
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
    Threads.@threads for j in 1:length(query_coord)
        ref_dist=[]
        for i in 1:length(ref_coord)
            ref_dist=push!(ref_dist, Distances.euclidean(query_coord[j],ref_coord[i]))
        end
        query_dist=push!(query_dist,findmin(ref_dist)[2])
        next!(p)
    end
    println("Polygons have been mapped to the cells!")
    my_annot=DataFrame(polygon_number=collect(1:length(query_dist)), mapped_cell=query_dist);
    from=sp.spmetaData.cell.cell
    to=sp.spmetaData.cell.cluster
    anno=map_values(my_annot, :mapped_cell, :cluster,from, to)
    sp.spmetaData.polygon=anno
    return sp
end
